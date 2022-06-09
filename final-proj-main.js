import { Articulated_Octopus } from "./arm.js";
import { tiny, defs } from "./examples/common.js";

// Pull these names into this module's scope for convenience:
const {
  vec3,
  vec4,
  color,
  Matrix,
  Mat4,
  Shape,
  Material,
  Shader,
  Texture,
  Component,
} = tiny;
import { Shape_From_File } from "./examples/obj-file-demo.js";

const Spline = class Spline {
  constructor() {
    this.points = [];
    this.tangents = [];
    this.size = 0;
  }

  add_point(x, y, z, tx, ty, tz) {
    this.points.push(vec3(x, y, z));
    this.tangents.push(vec3(tx, ty, tz));
    this.size = this.size + 1;
  }

  get_position(t) {
    if (this.size < 2) {
      return vec3(0, 0, 0);
    }
    const A = Math.floor(t * (this.size - 1));
    const B = Math.ceil(t * (this.size - 1));
    const s = (t * (this.size - 1)) % 1.0;

    let a = this.points[A].copy();
    let b = this.points[B].copy();
    return a.times(1 - s).plus(b.times(s));
  }

  set_point(index, x, y, z) {
    this.points[index] = vec3(x, y, z);
  }

  set_tangent(index, x, y, z) {
    this.tangents[index] = vec3(x, y, z);
  }

  reset() {
    this.points = [];
    this.tangents = [];
    this.size = 0;
  }
};

const Curved_Shape = class Curved_Shape extends Shape {
  constructor(curve_function, sample_count, curve_color = color(0, 1, 0, 1)) {
    super("position", "normal");

    this.material = {
      shader: new defs.Phong_Shader(),
      ambient: 1,
      color: curve_color,
    };
    this.sample_count = sample_count;

    this.curve_function = curve_function;
    if (curve_function && this.sample_count) {
      for (let i = 0; i < this.sample_count; i++) {
        let t = i / this.sample_count;
        this.arrays.position.push(curve_function(t));
        this.arrays.normal.push(vec3(0, 0, 0));
      }
    }
  }

  draw(webgl_manager, uniforms, transform) {
    super.draw(webgl_manager, uniforms, transform, this.material, "LINE_STRIP");
  }

  update(webgl_manager, uniforms, curve_function) {
    if (curve_function && this.sample_count) {
      for (let i = 0; i < this.sample_count; i++) {
        let t = i / this.sample_count;
        this.arrays.position[i] = curve_function(t);
      }
    }
  }

  get_arc_length() {
    let arcLen = 0;
    let prev = this.curve_function(0);
    for (let i = 1; i < this.sample_count; i++) {
      const t = i / this.sample_count;
      const curr = this.curve_function(t);
      arcLen = arcLen + curr.minus(prev).norm();
      prev = curr;
    }
    return arcLen;
  }
};

function h00(t) {
  return 2 * t ** 3 - 3 * t ** 2 + 1;
}

function h10(t) {
  return t ** 3 - 2 * t ** 2 + t;
}

function h01(t) {
  return -2 * t ** 3 + 3 * t ** 2;
}

function h11(t) {
  return t ** 3 - t ** 2;
}

function hermite(p0, p1, m0, m1, t) {
  const a = p0.times(h00(t));
  const b = m0.times(h10(t));
  const c = p1.times(h01(t));
  const d = m1.times(h11(t));
  return a.plus(b).plus(c).plus(d);
}

function compute_rotation(norm, diff_norm) {
  const v = norm.cross(diff_norm).normalized();
  const phi = -1 * Math.acos(norm.dot(diff_norm));
  const rcos = Math.cos(phi);
  const rsin = Math.sin(phi);
  const vx = v[0];
  const vy = v[1];
  const vz = v[2];
  const r = Matrix.of(
    [
      rcos + vx * vx * (1 - rcos),
      vz * rsin + vy * vx * (1 - rcos),
      -vy * rsin + vz * vx * (1 - rcos),
      0,
    ],
    [
      -vz * rsin + vx * vy * (1 - rcos),
      rcos + vy * vy * (1 - rcos),
      -vx * rsin + vz * vy * (1 - rcos),
      0,
    ],
    [
      vy * rsin + vx * vz * (1 - rcos),
      -vx * rsin + vy * vz * (1 - rcos),
      rcos + vz * vz * (1 - rcos),
      0,
    ],
    [0, 0, 0, 1]
  );
  return r;
}

const Hermite_Spline = class Hermite_Spline extends Spline {
  get_position(t) {
    if (this.size < 2) {
      return 0;
    }
    //do fancy shmancy stuff with h
    const A = Math.floor(t * (this.size - 1));
    const B = Math.ceil(t * (this.size - 1));
    const s = (t * (this.size - 1)) % 1.0;

    let tangentScale = 1 / (this.size - 1);
    let a = this.points[A].copy();
    let b = this.points[B].copy();
    let c = this.tangents[A].copy().times(tangentScale);
    let d = this.tangents[B].copy().times(tangentScale);
    return hermite(a, b, c, d, s);
  }
};

//x(t + dt) = x(t) + v(x(t), t) * dt
//v(t + dt) = v(t) + dt * (1/m) * f(x(t), v(t), t)
//Change the particle’s current position 𝐱􁈺𝑡􁈻 by
// integrating over the time step
function forward_euler(particle, dt, old_particle) {
  const mass = particle.mass;
  const force = particle.force;
  const prevVel = particle.velocity;
  const prevLoc = particle.position;
  //calculate v
  let newV = force.times(dt * (1 / mass)).plus(prevVel);
  particle.setVelocity(newV);
  //calculate x
  let newLoc = prevVel.times(dt).plus(prevLoc);
  particle.setPosition(newLoc);
}

//v(t + dt) = v(t) + dt * (1/m) * f(x(t), v(t), t)
//x(t + dt) = x(t) + dt * v(t + dt)
function sym_euler(particle, dt, old_particle) {
  const mass = particle.mass;
  const force = particle.force;
  const prevVel = particle.velocity;
  const prevLoc = particle.position;
  //let newV = prevVel + dt * (1/mass) * force;
  let newVel = force
    .times(1 / mass)
    .times(dt)
    .plus(prevVel);
  particle.setVelocity(newVel);
  let newLoc = prevLoc.plus(newVel.times(dt));
  particle.setPosition(newLoc);
}

//x(t + dt) = 2 * x(t) - x(t - dt) + (dt^2)/m * f(x(t), v(t), t)
//v(t) = (1/dt) * (x(t) - x(t-dt)
//velocity verlet: a = 1/m * f(x(t), v(t), t)
// x(t + dt) = x(t) + dt * v(t) + (dt^2)/2 * a(t)
// v(t + dt) = v(t) + dt/2 * (a(t) + a(t + dt))
function verlet_int(particle, dt, old_particle) {
  const mass = particle.mass;
  const prevVel = particle.velocity;
  const prevLoc = particle.position;
  const prevAccel = old_particle.force.times(1 / mass);
  const currAccel = particle.force.times(1 / mass);
  let newVel = prevAccel.plus(currAccel).times(0.5).times(dt).plus(prevVel);
  let newLoc = prevVel
    .times(dt)
    .plus(prevAccel.times(dt ** 2 / 2))
    .plus(prevLoc);
  particle.setVelocity(newVel);
  particle.setPosition(newLoc);
}

const Particle = class Particle {
  constructor(
    mass,
    position = vec3(0, 0, 0),
    spring_method,
    force = vec3(0, 0, 0),
    velocity = vec3(0, 0, 0)
  ) {
    this.mass = mass;
    this.position = position;
    this.force = force;
    this.velocity = velocity;
    this.spring_method = spring_method;
    this.transf = Mat4.identity();
  }

  zeroForce() {
    this.force = vec3(0, 0, 0);
  }

  apply_force(f) {
    this.force = this.force.plus(f);
  }

  setVelocity(v) {
    this.velocity = v;
  }

  setPosition(p) {
    this.position = p;
    this.transf = Mat4.translation(p[0], p[1], p[2]);
  }

  setMass(m) {
    this.mass = m;
  }
};

const Link = class Link {
  constructor() {
    //this indicates the indexes within the particle array that are linked by this object
    this.I = -1;
    this.J = -1;
    this.ks = -1;
    this.kd = -1;
    this.spring_length = -1;
    this.isEnabled = false;
  }

  link_particles(I, J, ks, kd, spring_length) {
    this.I = I;
    this.J = J;
    this.ks = ks;
    this.kd = kd;
    this.spring_length = spring_length;
    this.isEnabled = true;
  }
};

const Rock = class Rock {
  constructor(position, dimensions) {
    this.mass = 100;
    this.ks = 100;
    this.kd = 50;
    this.position = position;
    this.dimensions = dimensions;
  }
};

const Limb = class Limb {
  constructor(spring_method, transf, gravity, rocks) {
    this.transf = transf;
    this.particles = [];
    this.old_particles = [];
    this.links = [];
    this.spring_method = spring_method;
    this.rocks = rocks;

    //TODO make these varaibles configurable
    this.gravity = gravity;
    this.ground_ks = 3500.0;
    this.ground_kd = 10.0;
    this.dt = 0.03;
  }

  add_particle(mass, position) {
    this.particles.push(new Particle(mass, position, this.spring_method));
  }

  add_link(i, j, ks, kd, len) {
    this.links.push(new Link());
    this.links[this.links.length - 1].link_particles(i, j, ks, kd, len);
  }

  calculate_force(spring) {
    //assume this.particles[i] != this.particles[j] and both exist
    //assume spring is an initialized spring with ks, kd
    let j = spring.J;
    let i = spring.I;
    let len = spring.spring_length;
    //console.log(i + " " + j);
    const dij_vec = this.particles[j].position.minus(
      this.particles[i].position
    );
    const dij_mag =
      (dij_vec[0] ** 2 + dij_vec[1] ** 2 + dij_vec[2] ** 2) ** 0.5;
    const dij_hat = dij_vec.times(1 / dij_mag);
    const fs = dij_hat.times(dij_mag - len).times(spring.ks);
    const vij = this.particles[j].velocity.minus(this.particles[i].velocity);
    const fd = dij_hat.times(spring.kd).times(vij.dot(dij_hat));
    const fe = fs.plus(fd);
    const negFe = fe.times(-1);
    this.particles[j].apply_force(negFe);
    this.particles[i].apply_force(fe);
  }

  update(slows) {
    //copy all the particles before applying forces
    this.old_particles = [...this.particles];

    //zero out all particles' forces
    for (let i = 0; i < this.particles.length; i++) {
      this.particles[i].zeroForce();
    }

    //TODO: change to actual position of the limb
    //add position to the member variables and set particle position instead of topPos
    //const topPos = vec3(3, 0, 0);
    //this.particles[0].setPosition(topPos);

    //calculate spring movement and apply forces
    for (let i = 0; i < this.links.length; i++) {
      this.calculate_force(this.links[i]);
    }

    //apply force of gravity
    const gravForce = vec3(0, this.gravity, 0);
    for (let i = 1; i < this.particles.length; i++) {
      this.particles[i].apply_force(gravForce);
      //make water slow the speed a  lil -- fake friction
      if (slows) {
        this.particles[i].velocity = this.particles[i].velocity.times(0.99);
      }
    }

    for (let i = 0; i < this.particles.length; i++) {
      if (this.particles[i].position[1] < 0 && slows) {
        //console.log("i hit the ground owie");
        let len = 0.0;
        const dij_vec = vec3(0, -Math.abs(this.particles[i].position[1]), 0);
        const dij_mag =
          (dij_vec[0] ** 2 + dij_vec[1] ** 2 + dij_vec[2] ** 2) ** 0.5;
        const dij_hat = dij_vec.times(1 / dij_mag);
        const fs = dij_hat.times(dij_mag - len).times(this.ground_ks);
        const vij = vec3(0, 0, 0).minus(this.particles[i].velocity);
        const fd = dij_hat.times(this.ground_kd).times(vij.dot(dij_hat));
        const fe = fs.plus(fd.times(-1));
        const toApply = vec3(0, Math.abs(fe[1]), 0);
        //console.log(toApply);
        this.particles[i].apply_force(toApply);
        //friction -- stop after a ton of collisions
        this.particles[i].velocity = this.particles[i].velocity.times(0.99);
      }
    }
    for (let i = 0; i < this.rocks.length; i++) {
      this.handle_collisions(this.rocks[i]);
    }

    for (let i = 1; i < this.particles.length; i++) {
      this.spring_method(this.particles[i], this.dt, this.old_particles[i]);
    }
  }

  handle_collisions(rock) {
    for (let i = 0; i < this.particles.length; i++) {
      let pos = this.particles[i].position;
      let particlePos = this.transf.times(
        Mat4.translation(pos[0], pos[1], pos[2])
      );
      particlePos = vec3(
        particlePos[0][3],
        particlePos[1][3],
        particlePos[2][3]
      );
      if (
        particlePos[0] < rock.position[0] + rock.dimensions[0] &&
        particlePos[0] > rock.position[0] - rock.dimensions[0]
      ) {
        if (
          particlePos[1] < rock.position[1] + rock.dimensions[1] &&
          particlePos[1] > rock.position[1] - rock.dimensions[1]
        ) {
          if (
            particlePos[2] < rock.position[2] + rock.dimensions[2] &&
            particlePos[2] > rock.position[2] - rock.dimensions[2]
          ) {
            //we have collided! now determine which direction to hit the particle in...
            //console.log("collision!!");
            let diff = particlePos.minus(rock.position);
            diff = diff.times(10);
            this.particles[i].apply_force(diff);
          }
        }
      }
    }
  }
};

export const Final_Proj_base =
  (defs.Final_Proj_base = class Final_Proj_base extends Component {
    // **My_Demo_Base** is a Scene that can be added to any display canvas.
    // This particular scene is broken up into two pieces for easier understanding.
    // The piece here is the base class, which sets up the machinery to draw a simple
    // scene demonstrating a few concepts.  A subclass of it, Assignment2,
    // exposes only the display() method, which actually places and draws the shapes,
    // isolating that code so it can be experimented with on its own.
    init() {
      console.log("init");

      // constructor(): Scenes begin by populating initial values like the Shapes and Materials they'll need.
      this.hover = this.swarm = false;
      // At the beginning of our program, load one of each of these shape
      // definitions onto the GPU.  NOTE:  Only do this ONCE per shape it
      // would be redundant to tell it again.  You should just re-use the
      // one called "box" more than once in display() to draw multiple cubes.
      // Don't define more than one blueprint for the same thing here.
      this.shapes = {
        box: new defs.Cube(),
        ball: new defs.Subdivision_Sphere(4),
        axis: new defs.Axis_Arrows(),
        shell: new Shape_From_File("./assets/seashell.obj"),
        shell2: new Shape_From_File("./assets/shell2.obj"),
        octo: new Shape_From_File("./assets/octo.obj"),
        cave: new (defs.Subdivision_Sphere.prototype.make_flat_shaded_version())(2),
        starfish: new Shape_From_File("./assets/starfish.obj"),
        fish: new Shape_From_File("./assets/fish.obj"),
        cone : new defs.Cone_Tip ( 2, 10,  [[0.5,0],[0.5,0]] ),
        // cave: new defs.Subdivision_Sphere(2),
      };

      // *** Materials: ***  A "material" used on individual shapes specifies all fields
      // that a Shader queries to light/color it properly.  Here we use a Phong shader.
      // We can now tweak the scalar coefficients from the Phong lighting formulas.
      // Expected values can be found listed in Phong_Shader::update_GPU().
      const basic = new defs.Basic_Shader();
      const phong = new defs.Phong_Shader();
      const tex_phong = new defs.Textured_Phong();
      const my_water_shader = new defs.Water_Shader();
      this.materials = {};
      this.materials.plastic = {
        shader: phong,
        ambient: 0.2,
        diffusivity: 1,
        specularity: 0.5,
        color: color(0.9, 0.5, 0.9, 1),
      };
      this.materials.metal = {
        shader: phong,
        ambient: 0.2,
        diffusivity: 1,
        specularity: 1,
        color: color(0.9, 0.5, 0.9, 1),
      };
      this.materials.rgb = {
        shader: tex_phong,
        ambient: 0.5,
        texture: new Texture("assets/rgb.jpg"),
      };
      this.materials.water = {
        shader: my_water_shader,
        ambient: 0.2,
        diffusivity: 1,
        specularity: 1,
        color: color(0.1, 0.1, 0.9, 1),
      };

      this.materials.cave_texture = {
        shader: tex_phong,
        color: color(90/255, 90/255, 90/255, 1),
        ambient: 0.3, 
        diffusivity: 1, 
        specularity: .4,
        texture: new Texture("assets/rock.jpg")
      };

      this.materials.pink_coral = {
        shader: phong,
        ambient: 0.2,
        diffusivity: 1,
        specularity: 1,
        color: color(244/255, 194/255, 194/255, 1),
        texture: new Texture("assets/pinkcoral.jpeg")
      };

      this.materials.blue_coral = {
        shader: phong,
        ambient: 0.2,
        diffusivity: 1,
        specularity: 1,
        color: color(173/255, 216/255, 230/255, 1),
        texture: new Texture("assets/bluecoral.jpeg")
      };

      this.audio = new Audio("assets/somethingfishy.mp3");

      //Limb implementation
      this.spring_method = (p, t, x) => sym_euler(p, t, x);
      let ks = 8.9;
      let kd = 20.9;
      let len = 0.5;

      const seaweed_kd = 5.9;
      const seaweed_ks = 8.9;

      //assume meters
      this.gravity = -6.8;
      this.ground_ks = 350.0;
      this.ground_kd = 10.0;
      this.dt = 0.02;
      this.limbs = [];
      this.seaweed = [];
      this.isRunning = true;
      this.octopusSpeed = 0.15;

      //Rocks -- must be made FIRST
      this.rocks = [];
      this.rocks.push(new Rock(vec3(15, 0, -15), vec3(4.2, 4.5, 4.2)));
      this.rocks.push(new Rock(vec3(-16, 2, 0), vec3(2.8, 4, 4)));
      this.rocks.push(new Rock(vec3(0, -8, -30), vec3(50, 8, 50)));

      const LimbTransformArray = [
        Mat4.translation(0, -2, 3),
        Mat4.translation(3, -2, 0),
        Mat4.translation(-3, -2, 0),
        Mat4.translation(0, -2, -3),
        Mat4.translation(1.5, -2, 1.5),
        Mat4.translation(1.5, -2, -1.5),
        Mat4.translation(-1.5, -2, 1.5),
        Mat4.translation(-1.5, -2, -1.5),
      ];

      for (let i = 0; i < LimbTransformArray.length; i++) {
        this.limbs.push(
          new Limb(
            this.spring_method,
            LimbTransformArray[i],
            this.gravity * 0.5,
            this.rocks
          )
        );

        this.limbs[i].add_particle(1, vec3(0, 12, 0));
        this.limbs[i].add_particle(
          1,
          vec3(Math.random() * 8 - 4, 11, Math.random() * 8 - 4)
        );
        this.limbs[i].add_particle(
          1,
          vec3(Math.random() * 8 - 4, 10, Math.random() * 8 - 4)
        );
        this.limbs[i].add_particle(
          1,
          vec3(Math.random() * 8 - 4, 9, Math.random() * 8 - 4)
        );
        this.limbs[i].add_link(0, 1, ks * 2, kd * 0.5, len);
        this.limbs[i].add_link(1, 2, ks, kd, len * 1.5);
        this.limbs[i].add_link(2, 3, ks * 0.5, kd, len * 2.5);
      }
      //octopus motion setup
      //this.octopusPosition = Mat4.identity();
      this.octopusPosition = Mat4.translation(0, 14, 0);
      this.octopusDirection = Mat4.identity();

      const SeaweedTransform = [
        Mat4.translation(0, 0, 0),
        Mat4.translation(0, 0, 0),
        Mat4.translation(0, 0, 0),
        Mat4.translation(0, 0, 0),
        Mat4.translation(0, 0, 0),
      ];

      for (let i = 0; i < 5; i++) {
        let n = 13;
        this.seaweed.push(
          new Limb(
            this.spring_method,
            SeaweedTransform[i],
            this.gravity * -0.4,
            []
          )
        );
        this.seaweed[i].add_particle(1.5, vec3(0, 0, 0));
        this.seaweed[i].add_particle(
          0.5,
          vec3(Math.random() * 8 - 4, n, Math.random() * 8 - 4)
        );
        this.seaweed[i].add_particle(
          0.5,
          vec3(Math.random() * 8 - 4, n + 1, Math.random() * 8 - 4)
        );
        this.seaweed[i].add_particle(
          0.5,
          vec3(Math.random() * 8 - 4, n + 2, Math.random() * 8 - 4)
        );
        this.seaweed[i].add_particle(
          0.5,
          vec3(Math.random() * 8 - 4, n + 3, Math.random() * 8 - 4)
        );
        this.seaweed[i].add_link(0, 1, seaweed_ks, seaweed_kd, len / 2);
        this.seaweed[i].add_link(1, 2, seaweed_ks, seaweed_kd, len / 2);
        this.seaweed[i].add_link(2, 3, seaweed_ks, seaweed_kd, len / 2);
        this.seaweed[i].add_link(3, 4, seaweed_ks, seaweed_kd, len / 2);
      }

      this.spline = new Hermite_Spline();
      this.sample_count = 1000;
      const tangentScale = 10;
      let dx = 11;
      let dy = 3;
      let dz = -18;
      //Note: must be done to avoid WebGL complaints
      this.spline.add_point(0 + dx, 3 + dy, 4 + dz, 0, 0, 5);
      this.spline.add_point(4 + dx, 3 + dy, 8 + dz, 5, 0, 0);
      this.spline.add_point(8 + dx, 3 + dy, 4 + dz, 0, 0, -5);
      this.spline.add_point(4 + dx, 3 + dy, 0 + dz, -5, 0, 0);
      this.spline.add_point(0 + dx, 3 + dy, 4 + dz, 0, 0, 5);

      //const curve_func = (t) => vec3(0,0,0);
      const curve_func = (t) => this.spline.get_position(t);
      this.curve = new Curved_Shape(
        curve_func,
        this.sample_count,
        color(1, 0, 0, 1)
      );
      console.log(this.curve);

      // for ik octopus
      this.ik = false;
      this.octopus = null;
      this.limb_positions = [];
      this.limb_transf = [];
    }

    computeMovement(target) {
      if (this.octopus != null) {
        this.octopus.updateOctopus(target, 0);
      }
    }

    render_animation(caller) {
      // display():  Called once per frame of animation.  We'll isolate out
      // the code that actually draws things into Assignment2, a
      // subclass of this Scene.  Here, the base class's display only does
      // some initial setup.

      // Setup -- This part sets up the scene's overall camera matrix, projection matrix, and lights:
      if (!caller.controls) {
        this.animated_children.push(
          (caller.controls = new defs.Movement_Controls({
            uniforms: this.uniforms,
          }))
        );
        caller.controls.add_mouse_controls(caller.canvas);

        // Define the global camera and projection matrices, which are stored in shared_uniforms.  The camera
        // matrix follows the usual format for transforms, but with opposite values (cameras exist as
        // inverted matrices).  The projection matrix follows an unusual format and determines how depth is
        // treated when projecting 3D points onto a plane.  The Mat4 functions perspective() or
        // orthographic() automatically generate valid matrices for one.  The input arguments of
        // perspective() are field of view, aspect ratio, and distances to the near plane and far plane.

        // !!! Camera changed here
        // TODO: you can change the camera as needed.
        Shader.assign_camera(
          Mat4.look_at(vec3(5, 10, 30), vec3(0, 10, -10), vec3(0, 1, 0)),
          this.uniforms
        );
      }
      this.uniforms.projection_transform = Mat4.perspective(
        Math.PI / 4,
        caller.width / caller.height,
        1,
        100
      );

      // *** Lights: *** Values of vector or point lights.  They'll be consulted by
      // the shader when coloring shapes.  See Light's class definition for inputs.
      const t = (this.t = this.uniforms.animation_time / 1000);

      // const light_position = Mat4.rotation( angle,   1,0,0 ).times( vec4( 0,-1,1,0 ) ); !!!
      // !!! Light changed here
      const light_position = vec4(20, 20, 20, 1.0);
      this.uniforms.lights = [
        defs.Phong_Shader.light_source(
          light_position,
          color(1, 1, 1, 1),
          1000000
        ),
      ];

      // draw axis arrows.
      //this.shapes.axis.draw(caller, this.uniforms, Mat4.identity(), this.materials.rgb);
      this.audio.play();
    }
  });

export class Final_Proj extends Final_Proj_base {
  // **Assignment2** is a Scene object that can be added to any display canvas.
  // This particular scene is broken up into two pieces for easier understanding.
  // See the other piece, My_Demo_Base, if you need to see the setup code.
  // The piece here exposes only the display() method, which actually places and draws
  // the shapes.  We isolate that code so it can be experimented with on its own.
  // This gives you a very small code sandbox for editing a simple scene, and for
  // experimenting with matrix transformations.
  render_animation(caller) {
    // display():  Called once per frame of animation.  For each shape that you want to
    // appear onscreen, place a .draw() call for it inside.  Each time, pass in a
    // different matrix value to control where the shape appears.

    // Variables that are in scope for you to use:
    // this.shapes.box:   A vertex array object defining a 2x2x2 cube.
    // this.shapes.ball:  A vertex array object defining a 2x2x2 spherical surface.
    // this.materials.metal:    Selects a shader and draws with a shiny surface.
    // this.materials.plastic:  Selects a shader and draws a more matte surface.
    // this.lights:  A pre-made collection of Light objects.
    // this.hover:  A boolean variable that changes when the user presses a button.
    // shared_uniforms:  Information the shader needs for drawing.  Pass to draw().
    // caller:  Wraps the WebGL rendering context shown onscreen.  Pass to draw().

    // Call the setup code that we left inside the base class:
    super.render_animation(caller);

    /**********************************
         Start coding down here!!!!
         **********************************/
    // From here on down it's just some example shapes drawn for you -- freely
    // replace them with your own!  Notice the usage of the Mat4 functions
    // translation(), scale(), and rotation() to generate matrices, and the
    // function times(), which generates products of matrices.

    const blue = color(0, 0, 1, 1),
      yellow = color(1, 0.7, 0, 1),
      red = color(1, 0, 0, 1);
    const sand = color(211 / 255, 199 / 255, 162 / 255, 1);
    const ocean = color(173/255, 126 / 255, 230 / 255, .9);
    const lightShellColor = color(226 / 255, 223 / 255, 210 / 255, 0.75);
    const darkShellColor = color(247 / 255, 200 / 255, 194 / 255, 0.75);
    const starfishColor = color(250 / 255, 0 / 255, 127 / 255, 0.75);
    const seaweedColor = color(60 / 255, 130 / 255, 80 / 255, 1);
    const octoColor = color(135 / 255, 81 / 255, 109 / 255, 0.5);
    const white = color(1, 1, 1, 1);
    const black = color(0, 0, 0, 1);

    const t = (this.t = this.uniforms.animation_time / 1000);

    // !!! Draw ground
    let floor_transform = Mat4.translation(0, 0, 0).times(
      Mat4.scale(100, 0.01, 100)
    );
    this.shapes.box.draw(caller, this.uniforms, floor_transform, {
      ...this.materials.water,
    });

    //skybox
    let skybox_transform = Mat4.scale(50, 50, 50);
    this.shapes.ball.draw(caller, this.uniforms, skybox_transform, {
      ...this.materials.plastic,
      color: ocean,
    });

    //rocks
    for (let i = 0; i < this.rocks.length; i++) {
      let rock_transform = Mat4.translation(
        this.rocks[i].position[0],
        this.rocks[i].position[1],
        this.rocks[i].position[2]
      );
      rock_transform = rock_transform.times(
        Mat4.scale(
          this.rocks[i].dimensions[0],
          this.rocks[i].dimensions[1],
          this.rocks[i].dimensions[2]
        )
      );
      // this.shapes.box.draw(caller, this.uniforms, rock_transform, this.materials.cave_texture);
    }

    //random shells
    let firstShellTransform = Mat4.translation(-5, 1, -15)
      .times(Mat4.rotation(-Math.PI / 2, 0, 0, 1))
      .times(Mat4.rotation(-Math.PI / 2, 1, 0, 0))
      .times(Mat4.scale(.8, .8, .8));
    this.shapes.shell.draw(caller, this.uniforms, firstShellTransform, {
      ...this.materials.metal,
      color: lightShellColor,
    });

    let secondShellTransform = Mat4.translation(-15, 1, -20)
      .times(Mat4.rotation(-Math.PI / 2, 0, 0, 1))
      .times(Mat4.rotation(-Math.PI / 2, 1, 0, 0))
      .times(Mat4.scale(1.1, 1.1, 1.1));
    this.shapes.shell2.draw(caller, this.uniforms, secondShellTransform, {
      ...this.materials.metal,
      color: darkShellColor,
    });

    let thirdShellTransform = Mat4.translation(18, 1, -9)
      .times(Mat4.rotation(-Math.PI / 2, 0, 0, 1))
      .times(Mat4.rotation(2*Math.PI / 3, 0, 1, 0))
      .times(Mat4.scale(.3, .3, .3));
    this.shapes.shell.draw(caller, this.uniforms, thirdShellTransform, {
      ...this.materials.metal,
      color: darkShellColor,
    });

    //starfish!
    let starfishTransform = Mat4.translation(10, 1, 0)
      .times(Mat4.rotation(-Math.PI / 2, 0, 0, 1))
      .times(Mat4.rotation(-Math.PI / 2, 0, 1, 0))
      .times(Mat4.scale(1.2, 1.2, 1.2));
    this.shapes.starfish.draw(caller, this.uniforms, starfishTransform, {
      ...this.materials.metal,
      color: starfishColor,
    });


    // draw fish on rock, topPos is location of spline
    const spline_transform = Mat4.identity();
    //this.curve.draw(caller, this.uniforms, spline_transform);

    var currT = ((this.t + 5) % 5) / 5;
    const topPos = this.spline.get_position(currT);

    const fish_transform = Mat4.translation(...topPos)
    .times(Mat4.rotation(Math.PI / 2, 0, 0, 1))
    .times(Mat4.rotation(Math.PI / 2, 0, 1, 0));

    this.shapes.fish.draw(caller, this.uniforms, fish_transform, {
      ...this.materials.plastic,
      color: ocean,
    });

    if (!this.ik) {
      // draw particles and arms of octopus
      for (let i = 0; i < this.limbs.length; i++) {
        for (let j = 0; j < this.limbs[i].particles.length; j++) {
          let particleTransform = this.limbs[i].transf.times(
            this.limbs[i].particles[j].transf.times(Mat4.scale(0.3, 0.3, 0.3))
          );
          this.shapes.ball.draw(caller, this.uniforms, particleTransform, {
            ...this.materials.metal,
            color: octoColor,
          });
        }
      }
      // draw limbs
      for (let i = 0; i < this.limbs.length; i++) {
        for (let j = 0; j < this.limbs[i].links.length; j++) {
          const index = this.limbs[i].links[j].I;
          const index2 = this.limbs[i].links[j].J;
          const origPt = this.limbs[i].particles[index].position;
          const finalPt = this.limbs[i].particles[index2].position;
          const diff = finalPt.minus(origPt);
          const newPos3 = origPt.plus(diff.times(0.5));

          // reference to rotate link between particles:
          // https://stackoverflow.com/questions/52189123/calculate-rotation-matrix-to-transform-one-vector-to-another
          const norm = vec3(0, 1, 0);
          const diff_norm = diff.normalized();

          const transfLimb = this.limbs[i].transf
            .times(Mat4.translation(...newPos3))
            .times(compute_rotation(norm, diff_norm))
            .times(Mat4.scale(0.1 * (4 - j), diff.norm() / 2, 0.1 * (4 - j)));

          this.shapes.box.draw(caller, this.uniforms, transfLimb, {
            ...this.materials.metal,
            color: octoColor,
          });
        }
      }
      this.octopusPosition = this.octopusPosition.times(this.octopusDirection);

      //apply octopus movement to all the limbs too
      for (let i = 0; i < this.limbs.length; i++) {
        let newPos = this.limbs[i].particles[0].position.plus(
          vec3(
            this.octopusDirection[0][3],
            this.octopusDirection[1][3],
            this.octopusDirection[2][3]
          )
        );
        this.limbs[i].particles[0].setPosition(newPos);
      }
      //update with forces and stuff for next frame
      for (let i = 0; i < this.limbs.length; i++) {
        this.limbs[i].update(true);
      }
    } else {
      if (!this.octopus) {
        // get location of all the arms and points

        // 8x3 matrix of Mat4 position of joints
        this.limb_positions = [];
        for (let i = 0; i < this.limbs.length; i++) {
          let limb_i_positions = [];
          for (let j = 0; j < this.limbs[i].particles.length; j++) {
            let new_pos = this.limbs[i].particles[j].transf;
            if (limb_i_positions.length == 0) {
              for (let k = 0; k < 3; k++) {
                new_pos[k][3] -= this.octopusPosition[k][3];
              }
            } else {
              for (let k = 0; k < 3; k++) {
                new_pos[k][3] -= this.octopusPosition[k][3];
              }
              for (let l = 0; l < limb_i_positions.length; l++)
                for (let k = 0; k < 3; k++) {
                  new_pos[k][3] -= limb_i_positions[l][k][3];
                }
            }
            limb_i_positions.push(new_pos);
          }
          this.limb_positions.push(limb_i_positions);
        }

        // get transformation matrices for boxes of limbs
        this.limb_transf = [];
        for (let i = 0; i < this.limbs.length; i++) {
          let limb_i_transf = [];
          for (let j = 0; j < this.limbs[i].links.length; j++) {
            const index = this.limbs[i].links[j].I;
            const index2 = this.limbs[i].links[j].J;
            const origPt = this.limbs[i].particles[index].position;
            const finalPt = this.limbs[i].particles[index2].position;
            const diff = finalPt.minus(origPt);
            const newPos3 = diff.times(0.5);

            // reference to rotate link between particles:
            // https://stackoverflow.com/questions/52189123/calculate-rotation-matrix-to-transform-one-vector-to-another
            const norm = vec3(0, 1, 0);
            const diff_norm = diff.normalized();

            const transfLimb = this.limbs[i].transf
              .times(Mat4.translation(...newPos3))
              .times(compute_rotation(norm, diff_norm))
              .times(Mat4.scale(0.1 * (4 - j), diff.norm() / 2, 0.1 * (4 - j)));
            limb_i_transf.push(transfLimb);
          }
          this.limb_transf.push(limb_i_transf);
        }
        console.log(this.limb_positions[0].length);
        console.log(this.limb_transf[0].length);

        this.octopus = new Articulated_Octopus(
          this.octopusPosition,
          this.limb_positions,
          this.limb_transf,
          this.shapes.octo
        );
        this.octopus.root.location_matrix = this.octopusPosition;
        console.log(this.octopus);
        this.ik = true;
      }

      this.octopus.draw(caller, this.uniforms, {
        ...this.materials.metal,
        color: octoColor,
      });

      this.computeMovement(topPos);
    }

    // draw seaweed
    for (let i = 0; i < this.seaweed.length; i++) {
      for (let j = 0; j < this.seaweed[i].links.length; j++) {
        const index = this.seaweed[i].links[j].I;
        const index2 = this.seaweed[i].links[j].J;
        const origPt = this.seaweed[i].particles[index].position;
        const finalPt = this.seaweed[i].particles[index2].position;
        const diff = finalPt.minus(origPt);
        const newPos = origPt.plus(diff.times(0.33));
        const newPos2 = origPt.plus(diff.times(0.66));
        const transf = Mat4.identity().times(
          Mat4.translation(newPos[0] - 10, newPos[1], newPos[2]).times(
            Mat4.scale(0.2, 0.5, 0.2)
          )
        );
        const transf2 = Mat4.identity().times(
          Mat4.translation(newPos2[0] - 10, newPos2[1], newPos2[2]).times(
            Mat4.scale(0.2, 0.5, 0.2)
          )
        );

        this.shapes.box.draw(caller, this.uniforms, transf, {
          ...this.materials.metal,
          color: seaweedColor,
        });
        this.shapes.box.draw(caller, this.uniforms, transf2, {
          ...this.materials.metal,
          color: seaweedColor,
        });

        const transf3 = this.limbs[i].transf.times(
          Mat4.translation(newPos[0] + 8, newPos[1], newPos[2] + 5).times(
            Mat4.scale(0.2, 0.5, 0.2)
          )
        );
        const transf4 = this.limbs[i].transf.times(
          Mat4.translation(newPos2[0] + 8, newPos2[1], newPos2[2] + 5).times(
            Mat4.scale(0.2, 0.5, 0.2)
          )
        );

        this.shapes.box.draw(caller, this.uniforms, transf3, {
          ...this.materials.metal,
          color: seaweedColor,
        });
        this.shapes.box.draw(caller, this.uniforms, transf4, {
          ...this.materials.metal,
          color: seaweedColor,
        });

        const transf5 = this.limbs[i].transf.times(
          Mat4.translation(newPos[0], newPos[1], newPos[2] - 17).times(
            Mat4.scale(0.2, 0.5, 0.2)
          )
        );
        const transf6 = this.limbs[i].transf.times(
          Mat4.translation(newPos2[0], newPos2[1], newPos2[2] - 17).times(
            Mat4.scale(0.2, 0.5, 0.2)
          )
        );

        this.shapes.box.draw(caller, this.uniforms, transf5, {
          ...this.materials.metal,
          color: seaweedColor,
        });
        this.shapes.box.draw(caller, this.uniforms, transf6, {
          ...this.materials.metal,
          color: seaweedColor,
        });
      }

      this.seaweed[i].update();
    }

    // draw the body shape of the octopus
    // note that this is done regardless of ik mode or not, jus that it won't move
    let torsoTransform = this.octopusPosition.times(Mat4.scale(4.8, 4.2, 4.8));
    this.shapes.octo.draw(caller, this.uniforms, torsoTransform, {
      ...this.materials.plastic,
      color: octoColor,
    });

    // draw eyeballs
    let rightEyeTransform = this.octopusPosition.times(
      Mat4.translation(1.4, -1, 4).times(Mat4.scale(1, 1, 1))
    );
    this.shapes.ball.draw(caller, this.uniforms, rightEyeTransform, {
      ...this.materials.plastic,
      color: white,
    });
    let rightIrisTransform = this.octopusPosition.times(
      Mat4.translation(1.4, -1.1, 4.5).times(Mat4.scale(0.7, 0.7, 0.7))
    );
    this.shapes.ball.draw(caller, this.uniforms, rightIrisTransform, {
      ...this.materials.plastic,
      color: black,
    });
    let leftEyeTransform = this.octopusPosition.times(
      Mat4.translation(-1.4, -1, 4).times(Mat4.scale(1, 1, 1))
    );
    this.shapes.ball.draw(caller, this.uniforms, leftEyeTransform, {
      ...this.materials.plastic,
      color: white,
    });
    let leftIrisTransform = this.octopusPosition.times(
      Mat4.translation(-1.4, -1.1, 4.5).times(Mat4.scale(0.7, 0.7, 0.7))
    );
    this.shapes.ball.draw(caller, this.uniforms, leftIrisTransform, {
      ...this.materials.plastic,
      color: black,
    });

    //draw cave

    let model_transform_cave1 = Mat4.translation(15, 2, -14).times(Mat4.scale(5.5, 3.5, 5.5));
    this.shapes.cave.draw(caller, this.uniforms, model_transform_cave1, this.materials.cave_texture);

    let model_transform_small1 = Mat4.translation(12, 0, -10).times(Mat4.scale(1.5, 1, 1.5));
    this.shapes.cave.draw(caller, this.uniforms, model_transform_small1, this.materials.cave_texture);

    let model_transform_small2 = Mat4.translation(10.5, 0.4, -11).times(Mat4.scale(1.5, .75, 1));
    this.shapes.cave.draw(caller, this.uniforms, model_transform_small2, this.materials.cave_texture);

    let model_transform_small3 = Mat4.translation(20.5, 0.4, -13).times(Mat4.scale(2.5, 1.75, 2));
    this.shapes.cave.draw(caller, this.uniforms, model_transform_small3, this.materials.cave_texture);

    let model_transform_cave2 = Mat4.translation(-16, 2, 0).times(Mat4.scale(3.5, 4.5, 5.5));
    this.shapes.cave.draw(caller, this.uniforms, model_transform_cave2, this.materials.cave_texture);

    let model_transform_med1 = Mat4.translation(-12, 0, 0).times(Mat4.scale(1.75, 2, 1.5));
    this.shapes.cave.draw(caller, this.uniforms, model_transform_med1, this.materials.cave_texture);

    let model_transform_med2 = Mat4.translation(-15, 0.4, 4).times(Mat4.scale(1.5, 2, 1.5));
    this.shapes.cave.draw(caller, this.uniforms, model_transform_med2, this.materials.cave_texture);

    // draw coral
    let model_transform = Mat4.scale(5, 5, 5);

    let branch1_model = model_transform
      .times(Mat4.translation(-0.5, 0, 1))
      .times(Mat4.rotation(0.8, 0, 0, 1))
      .times(Mat4.scale(0.05, 0.4, 0.05));
    let branch2_model = branch1_model
      .times(Mat4.scale(10, 2, 10))
      .times(Mat4.translation(0, 1, 0))
      .times(Mat4.translation(-0.1, -0.5, 0))
      .times(Mat4.rotation(-0.7, 0, 0, 1))
      .times(Mat4.translation(0.1, 0.5, 0))
      .times(Mat4.scale(0.1, 0.5, 0.1));
    let branch3_model = model_transform
      .times(Mat4.translation(-0.4, -0.5, 1))
      .times(Mat4.rotation(0.1, 0, 0, 1))
      .times(Mat4.scale(0.05, 0.8, 0.05));
    let branch4_model = branch3_model
      .times(Mat4.scale(10, 1, 10))
      .times(Mat4.translation(0, 1, 0))
      .times(Mat4.translation(-0.05, -0.5, 0))
      .times(Mat4.rotation(0.5, 0, 0, 1))
      .times(Mat4.translation(0.25, 0.8, 0))
      .times(Mat4.scale(0.1, 0.5, 0.1));
    let branch5_model = branch3_model
      .times(Mat4.scale(10, 0.8, 10))
      .times(Mat4.translation(0, 1, 0))
      .times(Mat4.translation(-0.05, -0.5, 0))
      .times(Mat4.rotation(-0.3, 0, 0, 1))
      .times(Mat4.translation(-0.1, 1, 0))
      .times(Mat4.scale(0.1, 0.5, 0.1));
    let branch6_model = model_transform
      .times(Mat4.translation(-5.2, -2.4, 1))
      .times(Mat4.rotation(-0.2, 0, 0, 1))
      .times(Mat4.scale(0.05, 0.5, 0.05));
    let branch7_model = branch6_model
      .times(Mat4.scale(10, 2, 10))
      .times(Mat4.translation(0, 1, 0))
      .times(Mat4.translation(-0.05, -0.5, 0))
      .times(Mat4.rotation(-0.4, 0, 0, 1))
      .times(Mat4.translation(0.07, 0.5, 0))
      .times(Mat4.scale(0.1, 0.5, 0.1));
    let branch8_model = branch6_model
      .times(Mat4.scale(10, 1, 10))
      .times(Mat4.translation(0, 1, 0))
      .times(Mat4.translation(-0.05, -0.5, 0))
      .times(Mat4.rotation(0.8, 0, 0, 1))
      .times(Mat4.translation(1, 1.2, 0))
      .times(Mat4.scale(0.1, 0.5, 0.1));
    let branch9_model = model_transform
      .times(Mat4.translation(-5, -2.45, 1))
      .times(Mat4.rotation(-0.5, 0, 0, 1))
      .times(Mat4.scale(0.05, 0.5, 0.05));

    this.shapes.box.draw(
      caller,
      this.uniforms,
      branch1_model,
      this.materials.pink_coral
    );
    this.shapes.box.draw(
      caller,
      this.uniforms,
      branch2_model,
      this.materials.pink_coral
    );
    this.shapes.box.draw(
      caller,
      this.uniforms,
      branch3_model,
      this.materials.pink_coral
    );
    this.shapes.box.draw(
      caller,
      this.uniforms,
      branch4_model,
      this.materials.pink_coral
    );
    this.shapes.box.draw(
      caller,
      this.uniforms,
      branch5_model,
      this.materials.pink_coral
    );
    this.shapes.box.draw(
      caller,
      this.uniforms,
      branch6_model,
      this.materials.pink_coral
    );
    this.shapes.box.draw(
      caller,
      this.uniforms,
      branch7_model,
      this.materials.pink_coral
    );
    this.shapes.box.draw(
      caller,
      this.uniforms,
      branch8_model,
      this.materials.pink_coral
    );
    this.shapes.box.draw(
      caller,
      this.uniforms,
      branch9_model,
      this.materials.pink_coral
    );

    model_transform = Mat4.scale(5, 5, 5);

    branch1_model = model_transform
      .times(Mat4.translation(3.5, 0, 1))
      .times(Mat4.rotation(0.8, 0, 0, 1))
      .times(Mat4.scale(0.05, 0.4, 0.05));
    branch2_model = branch1_model
      .times(Mat4.scale(10, 2, 10))
      .times(Mat4.translation(0, 1, 0))
      .times(Mat4.translation(-0.1, -0.5, 0))
      .times(Mat4.rotation(-0.7, 0, 0, 1))
      .times(Mat4.translation(0.1, 0.5, 0))
      .times(Mat4.scale(0.1, 0.5, 0.1));
    branch3_model = model_transform
      .times(Mat4.translation(3.6, -0.5, 1))
      .times(Mat4.rotation(0.1, 0, 0, 1))
      .times(Mat4.scale(0.05, 0.8, 0.05));
    branch4_model = branch3_model
      .times(Mat4.scale(10, 1, 10))
      .times(Mat4.translation(0, 1, 0))
      .times(Mat4.translation(-0.05, -0.5, 0))
      .times(Mat4.rotation(0.5, 0, 0, 1))
      .times(Mat4.translation(0.25, 0.8, 0))
      .times(Mat4.scale(0.1, 0.5, 0.1));
    branch5_model = branch3_model
      .times(Mat4.scale(10, 0.8, 10))
      .times(Mat4.translation(0, 1, 0))
      .times(Mat4.translation(-0.05, -0.5, 0))
      .times(Mat4.rotation(-0.3, 0, 0, 1))
      .times(Mat4.translation(-0.1, 1, 0))
      .times(Mat4.scale(0.1, 0.5, 0.1));
    branch6_model = model_transform
      .times(Mat4.translation(-5.2, -2.4, 1))
      .times(Mat4.rotation(-0.2, 0, 0, 1))
      .times(Mat4.scale(0.05, 0.5, 0.05));
    branch7_model = branch6_model
      .times(Mat4.scale(10, 2, 10))
      .times(Mat4.translation(0, 1, 0))
      .times(Mat4.translation(-0.05, -0.5, 0))
      .times(Mat4.rotation(-0.4, 0, 0, 1))
      .times(Mat4.translation(0.07, 0.5, 0))
      .times(Mat4.scale(0.1, 0.5, 0.1));
    branch8_model = branch6_model
      .times(Mat4.scale(10, 1, 10))
      .times(Mat4.translation(0, 1, 0))
      .times(Mat4.translation(-0.05, -0.5, 0))
      .times(Mat4.rotation(0.8, 0, 0, 1))
      .times(Mat4.translation(1, 1.2, 0))
      .times(Mat4.scale(0.1, 0.5, 0.1));
    branch9_model = model_transform
      .times(Mat4.translation(-5, -2.45, 1))
      .times(Mat4.rotation(-0.5, 0, 0, 1))
      .times(Mat4.scale(0.05, 0.5, 0.05));

    this.shapes.box.draw(
      caller,
      this.uniforms,
      branch1_model,
      this.materials.blue_coral
    );
    this.shapes.box.draw(
      caller,
      this.uniforms,
      branch2_model,
      this.materials.blue_coral
    );
    this.shapes.box.draw(
      caller,
      this.uniforms,
      branch3_model,
      this.materials.blue_coral
    );
    this.shapes.box.draw(
      caller,
      this.uniforms,
      branch4_model,
      this.materials.blue_coral
    );
    this.shapes.box.draw(
      caller,
      this.uniforms,
      branch5_model,
      this.materials.blue_coral
    );
    this.shapes.box.draw(
      caller,
      this.uniforms,
      branch6_model,
      this.materials.blue_coral
    );
    this.shapes.box.draw(
      caller,
      this.uniforms,
      branch7_model,
      this.materials.blue_coral
    );
    this.shapes.box.draw(
      caller,
      this.uniforms,
      branch8_model,
      this.materials.blue_coral
    );
    this.shapes.box.draw(
      caller,
      this.uniforms,
      branch9_model,
      this.materials.blue_coral
    );
  }

  render_controls() {
    // render_controls(): Sets up a panel of interactive HTML elements, including
    // buttons with key bindings for affecting this scene, and live info readouts.
    this.control_panel.innerHTML += "Final Project: Octopus Motion";
    this.new_line();
    this.key_triggered_button("Debug", ["Shift", "D"], null);
    this.new_line();
    this.key_triggered_button(
      "Move octopus left",
      ["j"],
      () => {
        this.octopusDirection = Mat4.translation(-1 * this.octopusSpeed, 0, 0);
      },
      undefined,
      () => {
        this.octopusDirection = Mat4.identity();
      }
    );
    this.new_line();
    this.key_triggered_button(
      "Move octopus forward",
      ["i"],
      () => {
        this.octopusDirection = Mat4.translation(
          0.0,
          0,
          -1 * this.octopusSpeed
        );
      },
      undefined,
      () => {
        this.octopusDirection = Mat4.identity();
      }
    );
    this.new_line();
    this.key_triggered_button(
      "Move octopus right",
      ["l"],
      () => {
        this.octopusDirection = Mat4.translation(this.octopusSpeed, 0, 0);
      },
      undefined,
      () => {
        this.octopusDirection = Mat4.identity();
      }
    );
    this.new_line();
    this.key_triggered_button(
      "Move octopus backward",
      ["k"],
      () => {
        this.octopusDirection = Mat4.translation(0, 0, this.octopusSpeed);
      },
      undefined,
      () => {
        this.octopusDirection = Mat4.identity();
      }
    );
    this.new_line();
    this.key_triggered_button(
      "Move octopus down",
      ["n"],
      () => {
        this.octopusDirection = Mat4.translation(0, -1 * this.octopusSpeed, 0);
      },
      undefined,
      () => {
        this.octopusDirection = Mat4.identity();
      }
    );
    this.new_line();
    this.key_triggered_button(
      "Move octopus up",
      ["m"],
      () => {
        this.octopusDirection = Mat4.translation(0, this.octopusSpeed, 0);
      },
      undefined,
      () => {
        this.octopusDirection = Mat4.identity();
      }
    );
    this.new_line();
    this.key_triggered_button("Feeding mode", ["v"], this.feeding_mode);
  }

  // change the current octopus into a set of 3DOF arms
  feeding_mode() {
    if (this.ik) {
      // reset octopus
      this.octopus = null;
      this.limb_positions = [];
      this.limb_transf = [];
    }
    this.ik = !this.ik;
  }
}
