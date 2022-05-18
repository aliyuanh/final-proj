import {tiny, defs} from './examples/common.js';

// Pull these names into this module's scope for convenience:
const { vec3, vec4, color, Mat4, Shape, Material, Shader, Texture, Component } = tiny;

const shapes = {
    'sphere': new defs.Subdivision_Sphere( 5 ),
};

export
const Articulated_Human = 
class Articulated_Human {
    constructor() {
        const sphere_shape = shapes.sphere;

        // torso node
        const torso_transform = Mat4.scale(1, 2.5, 0.5);
        this.torso_node = new Node("torso", sphere_shape, torso_transform);
        // root->torso
        const root_location = Mat4.translation(-2, 5, 0);
        this.root = new Arc("root", null, this.torso_node, root_location);

        // head node
        let head_transform = Mat4.scale(.6, .6, .6);
        head_transform.pre_multiply(Mat4.translation(0, .6, 0));
        this.head_node = new Node("head", sphere_shape, head_transform);
        // torso->neck->head
        const neck_location = Mat4.translation(0, 2.5, 0);
        this.neck = new Arc("neck", this.torso_node, this.head_node, neck_location);
        this.torso_node.children_arcs.push(this.neck);

        // right upper arm node
        let ru_arm_transform = Mat4.scale(1.2, .2, .2);
        ru_arm_transform.pre_multiply(Mat4.translation(1.2, 0, 0));
        this.ru_arm_node = new Node("ru_arm", sphere_shape, ru_arm_transform);
        // torso->r_shoulder->ru_arm
        const r_shoulder_location = Mat4.translation(0.6, 2, 0);
        this.r_shoulder = new Arc("r_shoulder", this.torso_node, this.ru_arm_node, r_shoulder_location);
        this.torso_node.children_arcs.push(this.r_shoulder)

        // right lower arm node
        let rl_arm_transform = Mat4.scale(1, .2, .2);
        rl_arm_transform.pre_multiply(Mat4.translation(1, 0, 0));
        this.rl_arm_node = new Node("rl_arm", sphere_shape, rl_arm_transform);
        // ru_arm->r_elbow->rl_arm
        const r_elbow_location = Mat4.translation(2.4, 0, 0);
        this.r_elbow = new Arc("r_elbow", this.ru_arm_node, this.rl_arm_node, r_elbow_location);
        this.ru_arm_node.children_arcs.push(this.r_elbow)

        // right hand node
        let r_hand_transform = Mat4.scale(.4, .3, .2);
        r_hand_transform.pre_multiply(Mat4.translation(0.4, 0, 0));
        this.r_hand_node = new Node("r_hand", sphere_shape, r_hand_transform);
        // rl_arm->r_wrist->r_hand
        const r_wrist_location = Mat4.translation(2, 0, 0);
        this.r_wrist = new Arc("r_wrist", this.rl_arm_node, this.r_hand_node, r_wrist_location);
        this.rl_arm_node.children_arcs.push(this.r_wrist)
    }

    draw(webgl_manager, uniforms, material) {
        this.matrix_stack = [];
        this._rec_draw(this.root, Mat4.identity(), webgl_manager, uniforms, material);
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

    debug(arc=null) {
        if (arc === null)
            arc = this.root;

        if (arc !== this.root) {
            arc.articulation_matrix = arc.articulation_matrix.times(Mat4.rotation(0.02, 0, 0, 1));
        }

        const node = arc.child_node;
        for (const next_arc of node.children_arcs) {
            this.debug(next_arc);
        }
    }
}

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