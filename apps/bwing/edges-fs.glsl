#version 330

layout( location = 0 ) out vec4 fragColor;

noperspective in vec3 v_Position;

uniform sampler2D image;

uniform samplerBuffer colormap;
uniform samplerBuffer field;
uniform samplerBuffer texcoord;

uniform vec3 constant_color;
uniform int use_constant_color;
uniform int have_tessellation_shader;
uniform int u_lighting;
uniform int u_clip;
uniform vec3 u_clip_center;
uniform vec3 u_clip_normal;
uniform int u_field_mode;
uniform int u_picking = -1;

uniform mat4 u_ModelViewMatrix;
uniform mat4 u_NormalMatrix;

// TODO: make these uniforms
const int ncolor = 256;

uniform float u_umin;
uniform float u_umax;

uniform int u_edges;
uniform float u_alpha;

noperspective in vec3 altitude;
in vec3 v_Normal;
flat in int id;
noperspective in vec2 v_Barycentric;

void
get_color( float u , inout vec3 color, in int alpha ) {

  float umin = u_umin;
  float umax = u_umax;

  int indx = int(ncolor * (u - umin) / (umax - umin));

  if (indx < 0) indx = 0;
  if (indx > 255) indx = 255;

  color = (1 - alpha) * color + alpha * texelFetch(colormap, indx).xyz;
}

void
shading( in vec3 l , in vec3 n , in vec3 color , out vec3 color_out ) {

  float diffuse = max(0.0,dot(l,n));
  float phong = 128.0;
  float specular = pow(max(0.0,dot(-reflect(l,n),n)),phong);

  vec3 cd = color * diffuse;
  vec3 cs = vec3(0.1) * specular;
  vec3 ca = vec3(0.2);
  color_out = ca + cd + cs;
}

void main() {

  //vec3 color = vec3(0.9, 0.9, 0.0);
  vec3 color = vec3(0.4, 0.2, 0.8);

  float u = texelFetch(field, id).r;
  get_color(u, color, u_field_mode);

  vec3 color_out = color;

  fragColor = vec4(color_out, 1.0);//u_alpha);

}
