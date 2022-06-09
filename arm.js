import { tiny, defs } from "./examples/common.js";

// Pull these names into this module's scope for convenience:
const { vec3, vec4, color, Mat4, Shape, Material, Shader, Texture, Component } =
  tiny;
import { Matrix, Vector } from "./linalg-src/linalg.js";

const shapes = {
  sphere: new defs.Subdivision_Sphere(5),
  box: new defs.Cube(),
};

export const Articulated_Octopus = class Articulated_Octopus {
  constructor(torso_loc, particles_transform, limbs_transform, octopus_shape) {
    const sphere_shape = shapes.sphere;
    const box_shape = shapes.box;

    // torso node
    const torso_transform = Mat4.scale(4.8, 4.2, 4.8);
    this.torso_node = new Node("torso", octopus_shape, torso_transform);
    // root->torso
    const root_location = torso_loc;
    this.root = new Arc("root", null, this.torso_node, root_location);

    this.body_arcs = [];
    this.upper_arcs = [];
    this.middle_arcs = [];

    this.upper_arms = [];
    this.middle_arms = [];
    this.lower_arms = [];

    for (let i = 0; i < particles_transform.length; i++) {
      let u_arm_transform = limbs_transform[i][0];
      let m_arm_transform = limbs_transform[i][1];
      let l_arm_transform = limbs_transform[i][2];

      let body_loc = particles_transform[i][0];
      let upper_loc = particles_transform[i][1];
      let mid_loc = particles_transform[i][2];
      let lower_loc = particles_transform[i][3]

      let upper_node = new Node("upper " + i, box_shape, u_arm_transform);
      let body_arc = new Arc(
        "body " + i,
        this.torso_node,
        upper_node,
        body_loc
      );

      let middle_node = new Node(
        "middle " + i,
        box_shape,
        m_arm_transform
      );

      let upper_arc = new Arc("upper " + i, upper_node, middle_node, upper_loc);
      let lower_node = new Node("lower " + i, box_shape, l_arm_transform);
      let middle_arc = new Arc("middle " + i, middle_node, lower_node, mid_loc);
      
      this.body_arcs.push(body_arc);
      this.upper_arcs.push(upper_arc);
      this.middle_arcs.push(middle_arc);

      this.upper_arms.push(upper_node);
      this.middle_arms.push(middle_node);
      this.lower_arms.push(lower_node);

      this.torso_node.children_arcs.push(body_arc);
      upper_node.children_arcs.push(upper_arc);
      middle_node.children_arcs.push(middle_arc);
    }

    this.allowedError = 0.01;
  }

  draw(webgl_manager, uniforms, material) {
    this.matrix_stack = [];
    this._rec_draw(
      this.root,
      Mat4.identity(),
      webgl_manager,
      uniforms,
      material
    );
  }

  updateArm(i, target, depth, thetas) {
    let upperRot = Mat4.rotation(thetas[9 * i + 2], 0, 0, 1)
      .times(Mat4.rotation(thetas[9 * i + 1], 0, 1, 0))
      .times(Mat4.rotation(thetas[9 * i], 1, 0, 0));

    let middleRot = Mat4.rotation(thetas[9 * i + 5], 0, 0, 1)
      .times(Mat4.rotation(thetas[9 * i + 4], 0, 1, 0))
      .times(Mat4.rotation(thetas[9 * i + 3], 1, 0, 0));

    let lowerRot = Mat4.rotation(thetas[9 * i + 8], 0, 0, 1)
      .times(Mat4.rotation(thetas[9 * i + 7], 0, 1, 0))
      .times(Mat4.rotation(thetas[9 * i + 6], 1, 0, 0));

    // NOTE: commented out various lines so observe individual angle changes
    
    // this.body_arcs[i].articulation_matrix =
    //   this.body_arcs[i].articulation_matrix.times(upperRot);

    this.upper_arcs[i].articulation_matrix =
      this.upper_arcs[i].articulation_matrix.times(middleRot);

    // this.middle_arcs[i].articulation_matrix =
    //   this.middle_arcs[i].articulation_matrix.times(lowerRot);

    // calcalate error and keep going
    let endEffector = this.getEndEffector(this.root, Mat4.identity(), i);
    const endEffPos = vec3(
      endEffector[0][3],
      endEffector[1][3],
      endEffector[2][3]
    );
    let error = target.minus(endEffPos);
    error = error.norm();

    //max iterations of this will be 10
    //I do not want to burn your computer
    if (error > this.allowedError && depth < 10) {
      this.updateArm(i, target, depth + 1, thetas);
    }
  }

  updateOctopus(target, depth) {
    let iden = Vector.create(1);
    let newThetas = this.computeJacobian(target);
    // let thetas = [0, 0, 0, 0, 0, 0, 0]; //7 0's
    let thetas = new Array(8 * 3 * 3).fill(0);
    for (let i = 0; i < thetas.length; i++) {
      thetas[i] = newThetas[i]
    }
    // go through each of 8 arms

    // NOTE: change end of i to 1 so only one arm moves
    for (let i = 0; i < 1; i++) {
      this.updateArm(i, target, depth, thetas);
    }
  }

  computeJacobianForArm(target, i) {
    let endEffector = this.getEndEffector(this.root, Mat4.identity(), i);
    const endEffPos = vec3(
      endEffector[0][3],
      endEffector[1][3],
      endEffector[2][3]
    );

    let pos = this.root.location_matrix;
    pos = pos.times(this.body_arcs[i].location_matrix);

    let body_arc_position = vec3(pos[0][3], pos[1][3], pos[2][3]);
    let r = endEffPos.minus(body_arc_position);
    let row1 = vec3(1, 0, 0).cross(r);
    let row2 = vec3(0, 1, 0).cross(r);
    let row3 = vec3(0, 0, 1).cross(r);

    pos = pos.times(this.upper_arcs[i].location_matrix);
    let upper_arc_position = vec3(pos[0][3], pos[1][3], pos[2][3]);
    r = endEffPos.minus(upper_arc_position);
    let row4 = vec3(1, 0, 0).cross(r);
    let row5 = vec3(0, 1, 0).cross(r);
    let row6 = vec3(0, 0, 1).cross(r);

    pos = pos.times(this.middle_arcs[i].location_matrix);
    let middle_arc_position = vec3(pos[0][3], pos[1][3], pos[2][3]);
    r = endEffPos.minus(middle_arc_position);
    let row7 = vec3(1, 0, 0).cross(r);
    let row8 = vec3(0, 1, 0).cross(r);
    let row9 = vec3(0, 0, 1).cross(r);

    let J = Matrix.create(
      [row1[0], row1[1], row1[2]],
      [row2[0], row2[1], row2[2]],
      [row3[0], row3[1], row3[2]],
      [row4[0], row4[1], row4[2]],
      [row5[0], row5[1], row5[2]],
      [row6[0], row6[1], row6[2]],
      [row7[0], row7[1], row7[2]],
      [row8[0], row8[1], row8[2]],
      [row9[0], row9[1], row9[2]]
    ).transpose

    let JPlus = J.pseudoInverse();

    let error = target.minus(endEffPos);
    let k = 0.15;
    let dx = error.times(k);
    let dxArray = Matrix.create([dx[0], dx[1], dx[2]]).transpose;
    let dTheta = JPlus.mult(dxArray);

    return dTheta;
  }

  computeJacobian(target) {
    this.dddkmatrix_stack = [];
    let dthetas = [];
    for (let i = 0; i < 8; i++) {
      dthetas.push(this.computeJacobianForArm(target, i));
    }

    let d = []
    for (let i = 0; i < dthetas.length; i++) {
      for (let j = 0; j < 9; j++) {
        d.push(dthetas[i][j][0]);
      }
    }
    return d
  }

  getEndEffector(arc, matrix, i) {
    //this.matrix_stack = [];
    let loc = this.root.location_matrix;
    let myarc = this.root.articulation_matrix;
    matrix.post_multiply(loc.times(myarc));
    loc = this.body_arcs[i].location_matrix;
    myarc = this.body_arcs[i].articulation_matrix;
    matrix.post_multiply(loc.times(myarc));
    loc = this.upper_arcs[i].location_matrix;
    myarc = this.upper_arcs[i].articulation_matrix;
    matrix.post_multiply(loc.times(myarc));
    loc = this.middle_arcs[i].location_matrix;
    myarc = this.middle_arcs[i].articulation_matrix;
    matrix.post_multiply(loc.times(myarc));
    // //account for hand length
    // matrix.post_multiply(Mat4.translation(0.5, 0, 0));
    // //console.log(matrix);
    return matrix;
  }

  _rec_draw(arc, matrix, webgl_manager, uniforms, material) {
    if (arc !== null) {
      const L = arc.location_matrix;
      const A = arc.articulation_matrix;
      matrix.post_multiply(L.times(A));
      this.matrix_stack.push(matrix.copy());

      const node = arc.child_node;
      const T = node.transform_matrix;
      matrix.post_multiply(T);
      node.shape.draw(webgl_manager, uniforms, matrix, material);

      matrix = this.matrix_stack.pop();
      for (const next_arc of node.children_arcs) {
        this.matrix_stack.push(matrix.copy());
        this._rec_draw(next_arc, matrix, webgl_manager, uniforms, material);
        matrix = this.matrix_stack.pop();
      }
    }
  }

  debug(arc = null) {
    if (arc === null) arc = this.root;

    if (arc !== this.root) {
      arc.articulation_matrix = arc.articulation_matrix.times(
        Mat4.rotation(0.02, 0, 0, 1)
      );
    }

    const node = arc.child_node;
    for (const next_arc of node.children_arcs) {
      this.debug(next_arc);
    }
  }
};

class Node {
  constructor(name, shape, transform) {
    this.name = name;
    this.shape = shape;
    this.transform_matrix = transform;
    this.children_arcs = [];
  }
}

class Arc {
  constructor(name, parent, child, location) {
    this.name = name;
    this.parent_node = parent;
    this.child_node = child;
    this.location_matrix = location;
    this.articulation_matrix = Mat4.identity();
  }
}
