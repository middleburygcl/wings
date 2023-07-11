#version 410

layout( location = 0 ) out vec4 fragColor;

noperspective in vec3 v_Position;

uniform samplerBuffer field;
uniform samplerBuffer colormap;

uniform vec3 constant_color;
uniform int use_constant_color;
uniform int have_tessellation_shader;
uniform int u_lighting;
uniform int u_clip;
uniform vec3 u_clip_center;
uniform vec3 u_clip_normal;
uniform int u_field_mode;

uniform mat4 u_ModelViewMatrix;
uniform mat4 u_NormalMatrix;

// TODO: make these uniforms
const int ncolor = 256;

uniform float u_umin = 0;
uniform float u_umax = 1;

uniform int u_edges;
uniform float u_alpha;

noperspective in vec3 altitude;
flat in int id;
in vec3 v_Normal;

void
get_color( float u , inout vec3 color ) {

  float umin = u_umin;
  float umax = u_umax;

  int indx = int(ncolor*(u - umin)/(umax - umin));

  if (indx < 0) indx = 0;
  if (indx > 255) indx = 255;

  color = (1 - u_field_mode) * color + u_field_mode * texelFetch(colormap, indx).xyz;
}

void
shading( in vec3 l , in vec3 n , in vec3 color , out vec3 color_out ) {

  float diffuse = max(0.0,dot(l,n));
  float phong = 128.0;
  float specular = pow(max(0.0,dot(-reflect(l,n),n)),phong);

  vec3 cd = color * diffuse;
  vec3 cs = vec3(0.2) * specular;
  vec3 ca = vec3(0.2);
  color_out = ca + cd + cs;
}

void main() {

  float d = min(min(altitude[0],altitude[1]),altitude[2]);
  float intensity = u_edges * exp2(-0.5 * d * d);

  vec3 position = normalize(v_Position);
  vec3 normal = normalize(v_Normal);

  vec3 color = vec3(0.8,0.8,0.8);
  float u = texelFetch(field, id).r;
  get_color(u, color);
  vec3 color_shaded;
  shading( -position , normal , color , color_shaded );

  vec3 color_out = color * ( 1 - u_lighting ) + u_lighting * color_shaded;

  fragColor = intensity * vec4(0, 0, 0, 1.0) + (1.0 - intensity) * vec4(color_out, u_alpha);
}
