import {tiny} from '../tiny-graphics.js';
// Pull these names into this module's scope for convenience:
const {Vector, Vector3, vec, vec3, vec4, color, Matrix, Mat4, Shape, Shader, Component} = tiny;

const defs = {};

export {tiny, defs};

const Basic_Shader = defs.Basic_Shader =
  class Basic_Shader extends Shader {
      // Basic_Shader is nearly the simplest way to subclass Shader, which stores and manages a GPU program.
      update_GPU (context, gpu_addresses, uniforms, model_transform, material) {
          // update_GPU():  Define how to synchronize our JavaScript's variables to the GPU's:
          const [P, C, M] = [uniforms.projection_transform, uniforms.camera_inverse, model_transform],
                PCM       = P.times (C).times (M);
          context.uniformMatrix4fv (gpu_addresses.projection_camera_model_transform, false,
                                    Matrix.flatten_2D_to_1D (PCM.transposed ()));
      }
      shared_glsl_code () {           // ********* SHARED CODE, INCLUDED IN BOTH SHADERS *********
          return `precision mediump float;
                  varying vec4 VERTEX_COLOR;
      `;
      }
      vertex_glsl_code () {          // ********* VERTEX SHADER *********
          return this.shared_glsl_code () + `
        attribute vec4 color;
        attribute vec3 position;                            // Position is expressed in object coordinates.
        uniform mat4 projection_camera_model_transform;

        void main() { 
          gl_Position = projection_camera_model_transform * vec4( position, 1.0 );      // Move vertex to final space.
          VERTEX_COLOR = color;                                 // Use the hard-coded color of the vertex.
        }`;
      }
      fragment_glsl_code () {         // ********* FRAGMENT SHADER *********
          return this.shared_glsl_code () + `
        void main() {                                                   
          gl_FragColor = VERTEX_COLOR;    // Directly use per-vertex colors for interpolation.
        }`;
      }
  };


const Funny_Shader = defs.Funny_Shader =
  class Funny_Shader extends Shader {
      update_GPU (context, gpu_addresses, uniforms, model_transform, material) {
          const [P, C, M] = [uniforms.projection_transform, uniforms.camera_inverse, model_transform],
                PCM       = P.times (C).times (M);
          context.uniformMatrix4fv (gpu_addresses.projection_camera_model_transform, false,
                                    Matrix.flatten_2D_to_1D (PCM.transposed ()));
          context.uniform1f (gpu_addresses.animation_time, uniforms.animation_time / 1000);
      }
      shared_glsl_code () {
          return `precision mediump float;
                  varying vec2 f_tex_coord;
      `;
      }
      vertex_glsl_code () {
          return this.shared_glsl_code () + `
        attribute vec3 position;                            // Position is expressed in object coordinates.
        attribute vec2 texture_coord;
        uniform mat4 projection_camera_model_transform;

        void main() {
          gl_Position = projection_camera_model_transform * vec4( position, 1.0 );  // Move vertex to final space
          f_tex_coord = texture_coord;                 // Supply the original texture coords for interpolation.
        }`;
      }
      fragment_glsl_code () {
          return this.shared_glsl_code () + `
        uniform float animation_time;
        void main() { 
          float a = animation_time, u = f_tex_coord.x, v = f_tex_coord.y;
          
          // To color in all pixels, use an arbitrary math function based only on time and UV texture coordinates.
          gl_FragColor = vec4(                                    
            2.0 * u * sin(17.0 * u ) + 3.0 * v * sin(11.0 * v ) + 1.0 * sin(13.0 * a),
            3.0 * u * sin(18.0 * u ) + 4.0 * v * sin(12.0 * v ) + 2.0 * sin(14.0 * a),
            4.0 * u * sin(19.0 * u ) + 5.0 * v * sin(13.0 * v ) + 3.0 * sin(15.0 * a),
            5.0 * u * sin(20.0 * u ) + 6.0 * v * sin(14.0 * v ) + 4.0 * sin(16.0 * a));
        }`;
      }
  };


const Phong_Shader = defs.Phong_Shader =
  class Phong_Shader extends Shader {
      constructor (num_lights = 2) {
          super ();
          this.num_lights = num_lights;
      }
      shared_glsl_code () {          // ********* SHARED CODE, INCLUDED IN BOTH SHADERS *********
          return ` 
        precision mediump float;
        const int N_LIGHTS = ` + this.num_lights + `;
        uniform float ambient, diffusivity, specularity, smoothness;
        uniform vec4 light_positions_or_vectors[N_LIGHTS], light_colors[N_LIGHTS];
        uniform float light_attenuation_factors[N_LIGHTS];
        uniform vec4 shape_color;
        uniform vec3 squared_scale, camera_center;

        varying vec3 N, vertex_worldspace;
                                             // ***** PHONG SHADING HAPPENS HERE: *****
        vec3 phong_model_lights( vec3 N, vec3 vertex_worldspace ) {
            vec3 E = normalize( camera_center - vertex_worldspace );
            vec3 result = vec3( 0.0 );
            for(int i = 0; i < N_LIGHTS; i++) {
                vec3 surface_to_light_vector = light_positions_or_vectors[i].xyz -
                                               light_positions_or_vectors[i].w * vertex_worldspace;
                float distance_to_light = length( surface_to_light_vector );

                vec3 L = normalize( surface_to_light_vector );
                vec3 H = normalize( L + E );
                
                  // Compute diffuse and specular components of Phong Reflection Model.
                float diffuse  =      max( dot( N, L ), 0.0 );
                float specular = pow( max( dot( N, H ), 0.0 ), smoothness );     // Use Blinn's "halfway vector" method.
                float attenuation = 1.0 / (1.0 + light_attenuation_factors[i] * distance_to_light * distance_to_light );


                vec3 light_contribution = shape_color.xyz * light_colors[i].xyz * diffusivity * diffuse
                                                          + light_colors[i].xyz * specularity * specular;

                result += attenuation * light_contribution;
              }
            return result;
          } `;
      }
      vertex_glsl_code () {           // ********* VERTEX SHADER *********
          return this.shared_glsl_code () + `
        attribute vec3 position, normal;                            // Position is expressed in object coordinates.

        uniform mat4 model_transform;
        uniform mat4 projection_camera_model_transform;

        void main() {                                                                
            gl_Position = projection_camera_model_transform * vec4( position, 1.0 );     // Move vertex to final space.
                                            // The final normal vector in screen space.
            N = normalize( mat3( model_transform ) * normal / squared_scale);

            vertex_worldspace = ( model_transform * vec4( position, 1.0 ) ).xyz;
          } `;
      }
      fragment_glsl_code () {          // ********* FRAGMENT SHADER *********
          return this.shared_glsl_code () + `
        void main() {                          
                                           // Compute an initial (ambient) color:
            gl_FragColor = vec4( shape_color.xyz * ambient, shape_color.w );
                                           // Compute the final color with contributions from lights:
            gl_FragColor.xyz += phong_model_lights( normalize( N ), vertex_worldspace );
          } `;
      }
      static light_source (position, color, size) {
          return {position, color, attenuation: 1 / size};
      }
      send_material (gl, gpu, material) {
          gl.uniform4fv (gpu.shape_color, material.color);
          gl.uniform1f (gpu.ambient, material.ambient);
          gl.uniform1f (gpu.diffusivity, material.diffusivity);
          gl.uniform1f (gpu.specularity, material.specularity);
          gl.uniform1f (gpu.smoothness, material.smoothness);
      }
      send_uniforms (gl, gpu, uniforms, model_transform) {
          const O = vec4 (0, 0, 0, 1), camera_center = uniforms.camera_transform.times (O).to3 ();
          gl.uniform3fv (gpu.camera_center, camera_center);

          // Use the squared scale trick from "Eric's blog" instead of inverse transpose matrix:
          const squared_scale = model_transform.reduce (
            (acc, r) => { return acc.plus (vec4 (...r).times_pairwise (r)); }, vec4 (0, 0, 0, 0)).to3 ();
          gl.uniform3fv (gpu.squared_scale, squared_scale);

          // Send the current matrices to the shader as a single pre-computed final matrix, the product.
          const PCM = uniforms.projection_transform.times (uniforms.camera_inverse).times (model_transform);
          gl.uniformMatrix4fv (gpu.model_transform, false, Matrix.flatten_2D_to_1D (model_transform.transposed ()));
          gl.uniformMatrix4fv (gpu.projection_camera_model_transform, false,
                               Matrix.flatten_2D_to_1D (PCM.transposed ()));

          if ( !uniforms.lights || !uniforms.lights.length)
              return;         // Lights omitted, ambient only

          const light_positions_flattened = [], light_colors_flattened = [];
          for (var i = 0; i < 4 * uniforms.lights.length; i++) {
              light_positions_flattened.push (uniforms.lights[ Math.floor (i / 4) ].position[ i % 4 ]);
              light_colors_flattened.push (uniforms.lights[ Math.floor (i / 4) ].color[ i % 4 ]);
          }
          gl.uniform4fv (gpu.light_positions_or_vectors, light_positions_flattened);
          gl.uniform4fv (gpu.light_colors, light_colors_flattened);
          gl.uniform1fv (gpu.light_attenuation_factors, uniforms.lights.map (l => l.attenuation));
      }
      update_GPU (context, gpu_addresses, uniforms, model_transform, material) {
          const defaults    = {color: color (0, 0, 0, 1), ambient: 0, diffusivity: 1, specularity: 1, smoothness: 40};
          let full_material = Object.assign (defaults, material);

          this.send_material (context, gpu_addresses, full_material);
          this.send_uniforms (context, gpu_addresses, uniforms, model_transform);
      }
  };


const Textured_Phong = defs.Textured_Phong =
  class Textured_Phong extends Phong_Shader {
      vertex_glsl_code () {         // ********* VERTEX SHADER *********
          return this.shared_glsl_code () + `
        varying vec2 f_tex_coord;
        attribute vec3 position, normal;                            // Position is expressed in object coordinates.
        attribute vec2 texture_coord;

        uniform mat4 model_transform;
        uniform mat4 projection_camera_model_transform;

        void main() {
            gl_Position = projection_camera_model_transform * vec4( position, 1.0 );     // Move vertex to final space.
                                              // The final normal vector in screen space.
            N = normalize( mat3( model_transform ) * normal / squared_scale);

            vertex_worldspace = ( model_transform * vec4( position, 1.0 ) ).xyz;
                                              // Turn the per-vertex texture coordinate into an interpolated variable.
            f_tex_coord = texture_coord;
          } `;
      }
      fragment_glsl_code () {        // ********* FRAGMENT SHADER *********
          return this.shared_glsl_code () + `
        varying vec2 f_tex_coord;
        uniform sampler2D texture;

        void main() {
            vec4 tex_color = texture2D( texture, f_tex_coord );       // Sample texture image in the correct place.
            if( tex_color.w < .01 ) discard;
                                                                     // Compute an initial (ambient) color:
            gl_FragColor = vec4( ( tex_color.xyz + shape_color.xyz ) * ambient, shape_color.w * tex_color.w );
                                                                     // Compute the final color with contributions from lights:
            gl_FragColor.xyz += phong_model_lights( normalize( N ), vertex_worldspace );
          } `;
      }
      update_GPU (context, gpu_addresses, uniforms, model_transform, material) {
          super.update_GPU (context, gpu_addresses, uniforms, model_transform, material);

          if (material.texture && material.texture.ready) {
              // Select texture unit 0 for the fragment shader Sampler2D uniform called "texture":
              context.uniform1i (gpu_addresses.texture, 0);
              // For this draw, use the texture image from correct the GPU buffer:
              material.texture.activate (context, 0);
          }
      }
  };


const Fake_Bump_Map = defs.Fake_Bump_Map =
  class Fake_Bump_Map extends Textured_Phong {
      fragment_glsl_code () {                            // ********* FRAGMENT SHADER *********
          return this.shared_glsl_code () + `
        varying vec2 f_tex_coord;
        uniform sampler2D texture;

        void main()  {        
            vec4 tex_color = texture2D( texture, f_tex_coord );       // Sample texture image in the correct place.
            if( tex_color.w < .01 ) discard;
                            
            // This time, slightly disturb normals based on sampling the same image that was used for texturing.
            vec3 bumped_N  = N + tex_color.rgb - .5*vec3(1,1,1);
            gl_FragColor = vec4( ( tex_color.xyz + shape_color.xyz ) * ambient, shape_color.w * tex_color.w );
            gl_FragColor.xyz += phong_model_lights( normalize( bumped_N ), vertex_worldspace );
          } `;
      }
  };

// TODO: you should implement the required classes here or in another file.
const Water_Shader = defs.Water_Shader =
    class Water_Shader extends Shader{
        send_material( gl, gpu, material )
        {                                       // send_material(): Send the desired shape-wide material qualities to the
            // graphics card, where they will tweak the Phong lighting formula.
            gl.uniform4fv( gpu.color,    material.color       );
            gl.uniform1f ( gpu.ambient,        material.ambient     );
            gl.uniform1f ( gpu.diffusivity,    material.diffusivity );
            gl.uniform1f ( gpu.specularity,    material.specularity );
            gl.uniform1f ( gpu.smoothness,     material.smoothness  );
        }


        update_GPU( context, gpu_addresses, program_state, model_transform, material )
        {
            // TODO (#EC 2): Pass the same information to the shader as for EC part 1.  Additionally
            // pass material.color to the shader.
            const [ P, C, M ] = [ program_state.projection_transform, program_state.camera_inverse, model_transform ],
                PCM = P.times( C ).times( M );
            context.uniformMatrix4fv( gpu_addresses.projection_camera_model_transform, false, Matrix.flatten_2D_to_1D( PCM.transposed() ) );
            context.uniform1f ( gpu_addresses.time, program_state.animation_time / 1000 );
            //context.uniform1f(gpu_addresses.color, material.color);
            const defaults = { color: color( 0,0,0,1 ), ambient: 0, diffusivity: 1, specularity: 1, smoothness: 40 };
            material = Object.assign( {}, defaults, material );
            this.send_material ( context, gpu_addresses, material );
        }
        // TODO (#EC 2):  Complete the shaders, displacing the input sphere's vertices as
        // a fireball effect and coloring fragments according to displacement.

        shared_glsl_code()            // ********* SHARED CODE, INCLUDED IN BOTH SHADERS *********
        { return `
            precision highp float;
            varying vec2 vUv;
            varying vec3 vPosition;
            varying float colorTransfer;
            float brightness = 22.;
            vec2 uvScale = vec2(.06,.06);
            uniform float time;
            float xScale = 1.;
            float yScale = 1.;
            float speed = .01;
            vec3 dotColor = vec3(52./255.,18./255.,99./255.);
            vec3 baseColor = vec3(1., 1.,1.);
            varying vec2 f_tex_coord;
              varying float disp;
              uniform vec4 sun_color;
              varying float noise;
              vec3 mod289(vec3 x) {
        return x - floor(x * (1.0 / 289.0)) * 289.0;
      }
      vec4 mod289(vec4 x) {
        return x - floor(x * (1.0 / 289.0)) * 289.0;
      }
      vec4 permute(vec4 x) {
        return mod289(((x*34.0)+1.0)*x);
      }
      vec4 taylorInvSqrt(vec4 r) {
        return 1.79284291400159 - 0.85373472095314 * r;
      }
      vec3 fade(vec3 t) {
        return t*t*t*(t*(t*6.0-15.0)+10.0);
      }
      // Classic Perlin noise
      float cnoise(vec3 P) {
        vec3 Pi0 = floor(P); // Integer part for indexing
        vec3 Pi1 = Pi0 + vec3(1.0); // Integer part + 1
        Pi0 = mod289(Pi0);
        Pi1 = mod289(Pi1);
        vec3 Pf0 = fract(P); // Fractional part for interpolation
        vec3 Pf1 = Pf0 - vec3(1.0); // Fractional part - 1.0
        vec4 ix = vec4(Pi0.x, Pi1.x, Pi0.x, Pi1.x);
        vec4 iy = vec4(Pi0.yy, Pi1.yy);
        vec4 iz0 = Pi0.zzzz;
        vec4 iz1 = Pi1.zzzz;
        vec4 ixy = permute(permute(ix) + iy);
        vec4 ixy0 = permute(ixy + iz0);
        vec4 ixy1 = permute(ixy + iz1);
        vec4 gx0 = ixy0 * (1.0 / 7.0);
        vec4 gy0 = fract(floor(gx0) * (1.0 / 7.0)) - 0.5;
        gx0 = fract(gx0);
        vec4 gz0 = vec4(0.5) - abs(gx0) - abs(gy0);
        vec4 sz0 = step(gz0, vec4(0.0));
        gx0 -= sz0 * (step(0.0, gx0) - 0.5);
        gy0 -= sz0 * (step(0.0, gy0) - 0.5);
        vec4 gx1 = ixy1 * (1.0 / 7.0);
        vec4 gy1 = fract(floor(gx1) * (1.0 / 7.0)) - 0.5;
        gx1 = fract(gx1);
        vec4 gz1 = vec4(0.5) - abs(gx1) - abs(gy1);
        vec4 sz1 = step(gz1, vec4(0.0));
        gx1 -= sz1 * (step(0.0, gx1) - 0.5);
        gy1 -= sz1 * (step(0.0, gy1) - 0.5);
        vec3 g000 = vec3(gx0.x,gy0.x,gz0.x);
        vec3 g100 = vec3(gx0.y,gy0.y,gz0.y);
        vec3 g010 = vec3(gx0.z,gy0.z,gz0.z);
        vec3 g110 = vec3(gx0.w,gy0.w,gz0.w);
        vec3 g001 = vec3(gx1.x,gy1.x,gz1.x);
        vec3 g101 = vec3(gx1.y,gy1.y,gz1.y);
        vec3 g011 = vec3(gx1.z,gy1.z,gz1.z);
        vec3 g111 = vec3(gx1.w,gy1.w,gz1.w);
        vec4 norm0 = taylorInvSqrt(vec4(dot(g000, g000), dot(g010, g010), dot(g100, g100), dot(g110, g110)));
        g000 *= norm0.x;
        g010 *= norm0.y;
        g100 *= norm0.z;
        g110 *= norm0.w;
        vec4 norm1 = taylorInvSqrt(vec4(dot(g001, g001), dot(g011, g011), dot(g101, g101), dot(g011, g011)));
        g001 *= norm1.x;
        g011 *= norm1.y;
        g101 *= norm1.z;
        g111 *= norm1.w;
        float n000 = dot(g000, Pf0);
        float n100 = dot(g100, vec3(Pf1.x, Pf0.yz));
        float n010 = dot(g010, vec3(Pf0.x, Pf1.y, Pf0.z));
        float n110 = dot(g110, vec3(Pf1.xy, Pf0.z));
        float n001 = dot(g001, vec3(Pf0.xy, Pf1.z));
        float n101 = dot(g101, vec3(Pf1.x, Pf0.y, Pf1.z));
        float n011 = dot(g011, vec3(Pf0.x, Pf1.yz));
        float n111 = dot(g111, Pf1);
        vec3 fade_xyz = fade(Pf0);
        vec4 n_z = mix(vec4(n000, n100, n010, n110), vec4(n001, n101, n011, n111), fade_xyz.z);
        vec2 n_yz = mix(n_z.xy, n_z.zw, fade_xyz.y);
        float n_xyz = mix(n_yz.x, n_yz.y, fade_xyz.x);
        return 2.2 * n_xyz;
      }
      float surface3 ( vec3 coord ) {
          float frequency = 7.0;
          float n = 0.4;
          n -= 1.0    * abs( cnoise( coord * frequency ) );
          n -= 1.5    * abs( cnoise( coord * frequency * 4.0 ) );
          n -= 1.25   * abs( cnoise( coord * frequency * 4.0 ) );
          return clamp( n, -0.6, 1.0 );
      }
      #define PI 3.141592653589793238462643383279
      `;
        }
        vertex_glsl_code()           // ********* VERTEX SHADER *********
        {  return this.shared_glsl_code() + `
      attribute vec2 texture_coord;
      attribute vec3 position;
      uniform mat4 projection_camera_model_transform;
      float turbulence( vec3 p ) {
          float t = -0.5;
          for (float f = 1.0 ; f <= 10.0 ; f++ ){
              float power = pow( 2.0, f );
              t += abs( cnoise( vec3( power * p ) ) / power );
          }
          return t;
      }
       float fireSpeed = .5;
      float pulseHeight = 0.0;
      float displacementHeight = .5;
      float turbulenceDetail = .7;
      attribute vec3 normal;
      attribute vec2 uv;
      attribute vec2 uv2;
     // void main() {
       //   noise = -0.8 * turbulence( turbulenceDetail * normal + ( time * .2 ) );
         // float b = pulseHeight * cnoise(
           //   0.05 * position + vec3( 1.0 * time )
          //);
          //float displacement = ( 0.0 - displacementHeight ) * noise + b;
          //float disp = displacement*30.;
          //vec3 newPosition = position + normal * displacement;
          //gl_Position = projection_camera_model_transform * vec4( newPosition, 1.0 );
      //}
        void main(){
          //vPosition = position;
          vUv = vec2((pow(position.x,1.)), pow(position.z,1.))/2.;
          vec2 uvMax = ( 2.0 * asin( sin( 2.0 * PI * vUv ) ) ) / PI;
          float n = surface3(vec3(uvMax * uvScale, time * speed));
           vec3 s = vec3( clamp( n, 0.0, 1.0 ) ) * dotColor * brightness;
           float mag = sqrt(s.x*s.x+s.y*s.y+s.z*s.z);
          vec4 newPosition = vec4(position.x, position.y*mag/5. + position.y,position.z, 1.0);
          gl_Position = projection_camera_model_transform*newPosition;
          //vec2 positionBoi = vUv;
          //float offset = time*speed;
          //float a = sin(3.14 * xScale  * positionBoi.x + offset) + cos(3.14 * 50.0 * positionBoi.y);
	      //float b = sin(3.14 * xScale * 5.0 * positionBoi.x) * 1.0;
	      //float c = sin(3.14 * yScale * position.y + offset * 5.0) + cos(3.14 * 50.0 * positionBoi.x);
	      //float d = sin(3.14 * yScale * 4.0 * positionBoi.y + offset * 1.0);
	      //float color = a+b+c+d;
          //vec3 newPosition = vec3(position.x, position.y*color, position.z);
          //colorTransfer = color*position.y;
          //gl_Position = projection_camera_model_transform * vec4(newPosition,1.0);
        }
       `;
        }
        fragment_glsl_code()           // ********* FRAGMENT SHADER *********
        { return this.shared_glsl_code() + `
    void main(){
        //vec2 position = vUv;
	     //float offset = time * speed;
	     //float a = sin(3.14 * xScale  * position.x + offset) + cos(3.14 * 50.0 * position.y);
	     //float b = sin(3.14 * xScale * 5.0 * position.x) * 1.0 ;
          //float c = sin(3.14 * yScale * position.y + offset * 5.0) + cos(3.14 * 50.0 * position.x);
          //float d = sin(3.14 * yScale * 4.0 * position.y + offset * 1.0);
          //float color = a + b + c + d;
	     //float diffR = baseColor.r - dotColor.r;
	     //float diffG = baseColor.g  - dotColor.g;
	     //float diffB = baseColor.b - dotColor.b;
	     //gl_FragColor = vec4(dotColor + vec3(diffR,diffG,diffB)*(color-.9)/7., 1.0);
        //gl_FragColor = vec4(1.,1.,1.,1.);
        vec2 uvMax = ( 2.0 * asin( sin( 2.0 * PI * vUv ) ) ) / PI;
        float n = surface3(vec3(uvMax * uvScale, time * speed));
        vec3 s = vec3( clamp( n, 0.0, 1.0 ) ) * dotColor * brightness;
        float diffX = baseColor.r - dotColor.r;
        float diffY = baseColor.g - dotColor.g;
        float diffZ = baseColor.b - dotColor.b;
        gl_FragColor = vec4(dotColor+n*vec3(.2-diffX,.2-diffY,.2-diffY)*.5, 1. );
       // float ratio = .7 + .5*sin(2.*3.14159*colorTransfer/5.);
	     // gl_FragColor = vec4(dotColor(colorTransfer+.1)*2.,colorTransfer*20.+.2);
    }
        `;
        }


    }
