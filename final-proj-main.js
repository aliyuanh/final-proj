import {tiny, defs} from './examples/common.js';

// Pull these names into this module's scope for convenience:
const { vec3, vec4, color, Mat4, Shape, Material, Shader, Texture, Component } = tiny;

// TODO: you should implement the required classes here or in another file.

const Spline =
    class Spline {
      constructor(){
        this.points = [];
        this.tangents = [];
        this.size = 0;
      }

      add_point(x, y, z, tx, ty, tz){
        this.points.push(vec3(x, y, z));
        this.tangents.push(vec3(tx, ty, tz));
        this.size = this.size + 1;
      }

      get_position(t){
        if(this.size < 2){
          return vec3(0, 0, 0);
        }
        const A = Math.floor(t * (this.size-1));
        const B = Math.ceil(t * (this.size-1));
        const s = (t * (this.size-1) % 1.0);

        let a = this.points[A].copy();
        let b = this.points[B].copy();
        return a.times(1-s).plus(b.times(s));
      }

      set_point(index, x, y, z){
        this.points[index] = vec3(x,y,z);
      }
      set_tangent(index, x, y, z){
        this.tangents[index] = vec3(x,y,z);
      }

      reset(){
        this.points = [];
        this.tangents = [];
        this.size = 0;
      }
    };

const Curved_Shape =
    class Curved_Shape extends Shape {
      constructor(curve_function, sample_count, curve_color = color(0,1,0,1)){
        super("position", "normal");

        this.material = {shader: new defs.Phong_Shader(), ambient: 1, color: curve_color};
        this.sample_count = sample_count;

        this.curve_function = curve_function;
        if(curve_function && this.sample_count){
          for(let i = 0; i < this.sample_count; i++){
            let t = i / this.sample_count;
            this.arrays.position.push(curve_function(t));
            this.arrays.normal.push(vec3(0,0,0));
          }
        }
      }

      draw(webgl_manager, uniforms, transform){
        super.draw(webgl_manager, uniforms, transform, this.material, "LINE_STRIP");
      }

      update(webgl_manager, uniforms, curve_function) {
        if(curve_function && this.sample_count) {
          for (let i = 0; i < this.sample_count; i++) {
            let t = i / this.sample_count;
            this.arrays.position[i] = curve_function(t);
          }
        }
      }

      get_arc_length(){
        let arcLen = 0;
        let prev = this.curve_function(0);
        for(let i = 1; i < this.sample_count; i++){
          const t = i / this.sample_count;
          const curr = this.curve_function(t);
          arcLen = arcLen + curr.minus(prev).norm();
          prev = curr;
        }
        return arcLen;
      }

    };

function h00(t){
  return 2 * t ** 3 - 3 * t ** 2 + 1;
}
function h10(t){
  return t**3 - 2 * t ** 2 + t;
}
function h01(t){
  return -2 * t ** 3 + 3 * t ** 2;
}
function h11(t){
  return t**3 - t**2;
}

function hermite(p0, p1, m0, m1, t){
  const a = p0.times(h00(t));
  const b = m0.times(h10(t));
  const c = p1.times(h01(t));
  const d = m1.times(h11(t));
  return a.plus(b).plus(c).plus(d);
}

const Hermite_Spline =
    class Hermite_Spline extends Spline {
      get_position(t) {
        if(this.size < 2){
          return 0;
        }
        //do fancy shmancy stuff with h
        const A = Math.floor(t * (this.size-1));
        const B = Math.ceil(t * (this.size-1));
        const s = (t * (this.size-1) % 1.0);

        let tangentScale = 1/(this.size - 1);
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
function forward_euler(particle, dt, old_particle){
  const mass = particle.mass;
  const force = particle.force;
  const prevVel = particle.velocity;
  const prevLoc = particle.position;
  //calculate v
  let newV = force.times(dt * (1/mass)).plus(prevVel);
  particle.setVelocity(newV);
  //calculate x
  let newLoc = prevVel.times(dt).plus(prevLoc);
  particle.setPosition(newLoc);
}

//v(t + dt) = v(t) + dt * (1/m) * f(x(t), v(t), t)
//x(t + dt) = x(t) + dt * v(t + dt)
function sym_euler(particle, dt, old_particle){
  const mass = particle.mass;
  const force = particle.force;
  const prevVel = particle.velocity;
  const prevLoc = particle.position;
  //let newV = prevVel + dt * (1/mass) * force;
  let newVel = force.times(1/mass).times(dt).plus(prevVel);
  particle.setVelocity(newVel);
  let newLoc = prevLoc.plus(newVel.times(dt));
  particle.setPosition(newLoc);
}

//x(t + dt) = 2 * x(t) - x(t - dt) + (dt^2)/m * f(x(t), v(t), t)
//v(t) = (1/dt) * (x(t) - x(t-dt)
//velocity verlet: a = 1/m * f(x(t), v(t), t)
// x(t + dt) = x(t) + dt * v(t) + (dt^2)/2 * a(t)
// v(t + dt) = v(t) + dt/2 * (a(t) + a(t + dt))
function verlet_int(particle, dt, old_particle){
  const mass = particle.mass;
  const prevVel = particle.velocity;
  const prevLoc = particle.position;
  const prevAccel = old_particle.force.times(1/mass);
  const currAccel = particle.force.times(1/mass);
  let newVel = prevAccel.plus(currAccel).times(.5).times(dt).plus(prevVel);
  let newLoc = prevVel.times(dt).plus(prevAccel.times(dt ** 2 / 2)).plus(prevLoc);
  particle.setVelocity(newVel);
  particle.setPosition(newLoc);
}

const Particle =
    class Particle{
      constructor(mass, position = vec3(0,0,0), spring_method, force = vec3(0,0,0), velocity = vec3(0,0,0)){
        this.mass = mass;
        this.position = position;
        this.force = force;
        this.velocity = velocity;
        this.spring_method = spring_method;
        this.transf = Mat4.identity();
      }

      zeroForce(){
        this.force = vec3(0,0,0);
      }

      apply_force(f){
        this.force = this.force.plus(f);
      }

      setVelocity(v){
        this.velocity = v;
      }

      setPosition(p){
        this.position = p;
        this.transf = Mat4.translation(p[0], p[1], p[2]);
      }

      setMass(m){
        this.mass = m;
      }

    };

const Link =
    class Link{
      constructor(){
        //this indicates the indexes within the particle array that are linked by this object
        this.I = -1;
        this.J = -1;
        this.ks = -1;
        this.kd = -1;
        this.spring_length = -1;
        this.isEnabled = false;

      }
      link_particles(I, J, ks, kd, spring_length){
        this.I = I;
        this.J = J;
        this.ks = ks;
        this.kd = kd;
        this.spring_length = spring_length;
        this.isEnabled = true;
      }
    };




export
const Final_Proj_base = defs.Final_Proj_base =
    class Final_Proj_base extends Component
    {                                          
      // **My_Demo_Base** is a Scene that can be added to any display canvas.
      // This particular scene is broken up into two pieces for easier understanding.
      // The piece here is the base class, which sets up the machinery to draw a simple
      // scene demonstrating a few concepts.  A subclass of it, Assignment2,
      // exposes only the display() method, which actually places and draws the shapes,
      // isolating that code so it can be experimented with on its own.
      init()
      {
        console.log("init")

        // constructor(): Scenes begin by populating initial values like the Shapes and Materials they'll need.
        this.hover = this.swarm = false;
        // At the beginning of our program, load one of each of these shape
        // definitions onto the GPU.  NOTE:  Only do this ONCE per shape it
        // would be redundant to tell it again.  You should just re-use the
        // one called "box" more than once in display() to draw multiple cubes.
        // Don't define more than one blueprint for the same thing here.
        this.shapes = { 'box'  : new defs.Cube(),
          'ball' : new defs.Subdivision_Sphere( 4 ),
          'axis' : new defs.Axis_Arrows() };

        // *** Materials: ***  A "material" used on individual shapes specifies all fields
        // that a Shader queries to light/color it properly.  Here we use a Phong shader.
        // We can now tweak the scalar coefficients from the Phong lighting formulas.
        // Expected values can be found listed in Phong_Shader::update_GPU().
        const basic = new defs.Basic_Shader();
        const phong = new defs.Phong_Shader();
        const tex_phong = new defs.Textured_Phong();
        this.materials = {};
        this.materials.plastic = { shader: phong, ambient: .2, diffusivity: 1, specularity: .5, color: color( .9,.5,.9,1 ) }
        this.materials.metal   = { shader: phong, ambient: .2, diffusivity: 1, specularity:  1, color: color( .9,.5,.9,1 ) }
        this.materials.rgb = { shader: tex_phong, ambient: .5, texture: new Texture( "assets/rgb.jpg" ) }

        this.ball_location = vec3(1, 1, 1);
        this.ball_radius = 0.25;


        // TODO: you should create the necessary shapes
        this.particles = [];
        this.old_particles = [];
        this.links = [];
        this.isRunning = false;

        //assume meters
        this.gravity = -4.8;
        this.ground_ks = 350.0;
        this.ground_kd = 10.0;
        this.spring_method = (p, t, x) => verlet_int(p, t, x);
        this.dt = 0.02;
        this.isRunning = true;

        //10 particles
        this.particles.push(new Particle(1, vec3(2,12,0), this.spring_method));
        this.particles.push(new Particle(1, vec3(0,11,1), this.spring_method));
        this.particles.push(new Particle(1, vec3(4,10,0), this.spring_method));
        this.particles.push(new Particle(1, vec3(0,9,-1), this.spring_method));

        this.particles.push(new Particle(1, vec3(-2,8,0), this.spring_method));
        this.particles.push(new Particle(1, vec3(0,7,2), this.spring_method));
        this.particles.push(new Particle(1, vec3(1,6,0), this.spring_method));
        this.particles.push(new Particle(1, vec3(0,5,2), this.spring_method));

        this.particles.push(new Particle(1, vec3(1,4,1), this.spring_method));
        this.particles.push(new Particle(1, vec3(0,3,0), this.spring_method));

        //9 links
        this.links.push(new Link());
        this.links.push(new Link());
        this.links.push(new Link());

        this.links.push(new Link());
        this.links.push(new Link());
        this.links.push(new Link());

        this.links.push(new Link());
        this.links.push(new Link());
        this.links.push(new Link());

        let ks = 8.9;
        let kd = 5.9;
        let len = .5;
        //this.gravity = -2.0;

        this.links[0].link_particles(0, 1, ks, kd, len);
        this.links[1].link_particles(1, 2, ks, kd, len);
        this.links[2].link_particles(2, 3, ks, kd, len);

        this.links[3].link_particles(3, 4, ks, kd, len);
        this.links[4].link_particles(4, 5, ks, kd, len);
        this.links[5].link_particles(5, 6, ks, kd, len);

        this.links[6].link_particles(6, 7, ks, kd, len);
        this.links[7].link_particles(7, 8, ks, kd, len);
        this.links[8].link_particles(8, 9, ks, kd, len);

        this.spline = new Hermite_Spline();
        this.sample_count = 1000;

        this.spline.add_point(3, 5, 0, -10, 0, 20);
        this.spline.add_point(1, 5, -2, 20, 0, 10);
        this.spline.add_point(2, 8, 5, 10, 0, -20);
        this.spline.add_point(-4, 5, 0, -20, 0, -10);
        this.spline.add_point(5, 6, -2, -10, 0, 20);

        const curve_func = (t) => this.spline.get_position(t);
        this.curve = new Curved_Shape(curve_func, this.sample_count);

      }

      render_animation( caller )
      {                                                // display():  Called once per frame of animation.  We'll isolate out
        // the code that actually draws things into Assignment2, a
        // subclass of this Scene.  Here, the base class's display only does
        // some initial setup.

        // Setup -- This part sets up the scene's overall camera matrix, projection matrix, and lights:
        if( !caller.controls )
        { this.animated_children.push( caller.controls = new defs.Movement_Controls( { uniforms: this.uniforms } ) );
          caller.controls.add_mouse_controls( caller.canvas );

          // Define the global camera and projection matrices, which are stored in shared_uniforms.  The camera
          // matrix follows the usual format for transforms, but with opposite values (cameras exist as
          // inverted matrices).  The projection matrix follows an unusual format and determines how depth is
          // treated when projecting 3D points onto a plane.  The Mat4 functions perspective() or
          // orthographic() automatically generate valid matrices for one.  The input arguments of
          // perspective() are field of view, aspect ratio, and distances to the near plane and far plane.

          // !!! Camera changed here
          // TODO: you can change the camera as needed.
          Shader.assign_camera( Mat4.look_at (vec3 (5, 8, 15), vec3 (0, 5, 0), vec3 (0, 1, 0)), this.uniforms );
        }
        this.uniforms.projection_transform = Mat4.perspective( Math.PI/4, caller.width/caller.height, 1, 100 );

        // *** Lights: *** Values of vector or point lights.  They'll be consulted by
        // the shader when coloring shapes.  See Light's class definition for inputs.
        const t = this.t = this.uniforms.animation_time/1000;

        // const light_position = Mat4.rotation( angle,   1,0,0 ).times( vec4( 0,-1,1,0 ) ); !!!
        // !!! Light changed here
        const light_position = vec4(20, 20, 20, 1.0);
        this.uniforms.lights = [ defs.Phong_Shader.light_source( light_position, color( 1,1,1,1 ), 1000000 ) ];

        // draw axis arrows.
        this.shapes.axis.draw(caller, this.uniforms, Mat4.identity(), this.materials.rgb);

        let spline_transform = Mat4.identity();
        this.curve.draw(caller, this.uniforms, spline_transform);
      }
    }


export class Final_Proj extends Final_Proj_base
{                                                    
  // **Assignment2** is a Scene object that can be added to any display canvas.
  // This particular scene is broken up into two pieces for easier understanding.
  // See the other piece, My_Demo_Base, if you need to see the setup code.
  // The piece here exposes only the display() method, which actually places and draws
  // the shapes.  We isolate that code so it can be experimented with on its own.
  // This gives you a very small code sandbox for editing a simple scene, and for
  // experimenting with matrix transformations.
  render_animation( caller )
  {                                                // display():  Called once per frame of animation.  For each shape that you want to
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
    super.render_animation( caller );

    /**********************************
     Start coding down here!!!!
     **********************************/
    // From here on down it's just some example shapes drawn for you -- freely
    // replace them with your own!  Notice the usage of the Mat4 functions
    // translation(), scale(), and rotation() to generate matrices, and the
    // function times(), which generates products of matrices.

    const blue = color( 0,0,1,1 ), yellow = color( 1,0.7,0,1 ), 
          wall_color = color( 0.7, 1.0, 0.8, 1 ), 
          blackboard_color = color( 0.2, 0.2, 0.2, 1 );
    const red = color(1, 0, 0, 1);

    const t = this.t = this.uniforms.animation_time/1000;

    // !!! Draw ground
    let floor_transform = Mat4.translation(0, 0, 0).times(Mat4.scale(10, 0.01, 10));
    this.shapes.box.draw( caller, this.uniforms, floor_transform, { ...this.materials.plastic, color: yellow } );

    // TODO: you should draw scene here.
    //let wall_transform = Mat4.translation(0, 5, -1.2).times(Mat4.scale(6, 5, 0.1));
    //this.shapes.box.draw( caller, this.uniforms, wall_transform, { ...this.materials.plastic, color: wall_color } );
    //let board_transform = Mat4.translation(3, 6, -1).times(Mat4.scale(2.5, 2.5, 0.1));
    //this.shapes.box.draw( caller, this.uniforms, board_transform, { ...this.materials.plastic, color: blackboard_color } );


    //all particles
    for(let i = 0; i < this.particles.length; i++){
      let particle_transform = this.particles[i].transf.times(Mat4.scale(.3, .3, .3));
      this.shapes.ball.draw(caller, this.uniforms, particle_transform, { ...this.materials.metal, color: red });
    }

    //copy all the particles before applying forces
    this.old_particles = [...this.particles];

    //zero out all particles' forces
    if(this.isRunning){
      for(let i = 0; i < this.particles.length; i++){
        this.particles[i].zeroForce();
      }
    }

    //do not apply external force to the top particle
    //calculate new position (sinusoidal based on time)
    //sin + 1 to shift above
    //then other particles should be able to calculate based on just the top particle moving??
    var currT = (Math.sin(this.t/3) + 1) / 2;
    currT = Math.min(Math.max(currT, 0.00001), .999999);
    const topPos = this.spline.get_position(currT);
    this.particles[0].setPosition(topPos);

    //calculate spring movement and apply forces
    for(let i = 0; i < this.links.length; i++){
      if(this.isRunning) {
        this.calculate_force(this.links[i]);
      }
      //draw balls between the links to act as spring visual indicators (.2 radius?)
      const index = this.links[i].I;
      const index2 = this.links[i].J;
      const origPt = this.particles[index].position;
      const finalPt = this.particles[index2].position;
      const diff = finalPt.minus(origPt);
      //draw 2 balls
      const newPos = origPt.plus(diff.times(0.33));
      const newPos2 = origPt.plus(diff.times(0.66));
      const transf = Mat4.translation(newPos[0], newPos[1], newPos[2]).times(Mat4.scale(.1, .1, .1));
      const transf2 = Mat4.translation(newPos2[0], newPos2[1], newPos2[2]).times(Mat4.scale(.1, .1, .1));

      this.shapes.ball.draw(caller, this.uniforms, transf, { ...this.materials.metal, color: blue });
      this.shapes.ball.draw(caller, this.uniforms, transf2, { ...this.materials.metal, color: blue });

    }

    //apply force of gravity
    const gravForce = vec3(0, this.gravity, 0);
    if(this.isRunning){
      for(let i = 1; i < this.particles.length; i++){
        this.particles[i].apply_force(gravForce);
      }
    }

    //check for collisions
    if(this.isRunning){
      for(let i = 1; i < this.particles.length; i++){
        if(this.particles[i].position[1] < 0){
          //going thru the ground
          let len = 0.0;
          const dij_vec = vec3(0, -Math.abs(this.particles[i].position[1]), 0);
          const dij_mag = (dij_vec[0] ** 2 + dij_vec[1] ** 2 + dij_vec[2] ** 2) ** 0.5;
          const dij_hat = dij_vec.times(1/dij_mag);
          const fs = dij_hat.times(dij_mag - len).times(this.ground_ks);
          const vij = vec3(0,0,0).minus(this.particles[i].velocity);
          const fd = dij_hat.times(this.ground_kd).times(vij.dot(dij_hat));
          const fe = fs.plus(fd.times(-1));
          const toApply = vec3(0, Math.abs(fe[1]), 0);
          this.particles[i].apply_force(toApply);
          //friction -- stop after a ton of collisions
          this.particles[i].velocity = this.particles[i].velocity.times(.99);

        }
      }
    }
    //calculate next positions using chosen movement type

    if(this.isRunning){
      for(let i = 1; i < this.particles.length; i++){
        this.spring_method(this.particles[i], this.dt, this.old_particles[i]);
      }
    }
  }

  render_controls()
  {                                 
    // render_controls(): Sets up a panel of interactive HTML elements, including
    // buttons with key bindings for affecting this scene, and live info readouts.
    this.control_panel.innerHTML += "Assignment 2: IK Engine";
    this.new_line();    
    // TODO: You can add your button events for debugging. (optional)
    this.key_triggered_button( "Debug", [ "Shift", "D" ], null );
    this.new_line();
  }

  //for every spring, add a force to the corresponding particles
  calculate_force(spring){
    //assume this.particles[i] != this.particles[j] and both exist
    //assume spring is an initialized spring with ks, kd
    let j = spring.J;
    let i = spring.I;
    let len = spring.spring_length;
    //console.log(i + " " + j);
    const dij_vec = this.particles[j].position.minus(this.particles[i].position);
    const dij_mag = (dij_vec[0] ** 2 + dij_vec[1] ** 2 + dij_vec[2] ** 2) ** 0.5;
    const dij_hat = dij_vec.times(1/dij_mag);
    const fs = dij_hat.times(dij_mag - len).times(spring.ks);
    const vij = this.particles[j].velocity.minus(this.particles[i].velocity);
    const fd = dij_hat.times(spring.kd).times(vij.dot(dij_hat));
    const fe = fs.plus(fd);
    const negFe = fe.times(-1);
    this.particles[j].apply_force(negFe);
    this.particles[i].apply_force(fe);
  }
}
