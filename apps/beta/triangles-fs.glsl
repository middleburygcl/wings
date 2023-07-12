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

uniform float u_umin = 0;
uniform float u_umax = 1;

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

  float diffuse = max(0.0,abs(dot(l,n)));
  float phong = 128.0;
  float specular = pow(max(0.0,dot(-reflect(l,n),n)),phong);

  vec3 cd = color * diffuse;
  vec3 cs = vec3(0.1) * specular;
  vec3 ca = vec3(0.2);
  color_out = ca + cd + cs;
}

void main() {

  float alpha = u_alpha;
  // if (u_alpha < 1.0) {
  //   float z = 1.0 / gl_FragCoord.w;
  //   alpha = u_alpha;//pow(u_alpha, 1./z);
  //   if (z < .1) {
  //     discard;
  //     return;
  //   }
  // }

  float d = min(min(altitude[0],altitude[1]),altitude[2]);
  float intensity = u_edges * exp2(-0.25 * d * d);

  vec3 position = normalize(v_Position);
  vec3 normal = normalize(v_Normal);

  vec3 color = vec3(0.8, 0.8, 0.8);

  #if ORDER == 0

  float u = texelFetch(field, id).r;
  get_color(u, color, u_field_mode);

  #elif ORDER == 1

  float u0 =  texelFetch(field, 3 * id    ).r;
  float u1 =  texelFetch(field, 3 * id + 1).r;
  float u2 =  texelFetch(field, 3 * id + 2).r;

  float s = v_Barycentric[0];
  float t = v_Barycentric[1];

  float u = (1.0 - s - t) * u0 + s * u1 + t * u2;
  get_color(u, color, u_field_mode);

  #elif ORDER == -1

  float s = v_Barycentric[0];
  float t = v_Barycentric[1];

  vec2 uv0 = texelFetch(texcoord, 3 * id    ).xy;
  vec2 uv1 = texelFetch(texcoord, 3 * id + 1).xy;
  vec2 uv2 = texelFetch(texcoord, 3 * id + 2).xy;

  vec2 uv = (1.0 - s - t) * uv0 + s * uv1 + t * uv2;

  float f = 200.0;
  if (sin(f * uv[0]) * sin(f * uv[1]) > 0.0)
    color = vec3(0, 0, 0);
  else
    color = vec3(1, 1, 1);

  // if we want parameter coordinates (checkerboard), u_field_mode = 1
  // if we want to plot the image, u_field_mode = 0
  color = (1 - u_field_mode) * texture(image, uv).rgb + u_field_mode * color;
  #endif

  vec3 color_shaded;
  shading( -position , normal , color , color_shaded );

  vec3 color_out = color * ( 1 - u_lighting ) + u_lighting * color_shaded;

  fragColor = intensity * vec4(0, 0, 0, 1.0) + (1.0 - intensity) * vec4(color_out, alpha);

  // modify the color if this element was picked
  // u_picking = -1 if not picking, u_picking = element id if picking
  if (u_picking == id) fragColor = intensity * vec4(0, 0, 0, 1.0) + (1.0 - intensity) * vec4(0.2, 0.5, 1.0, 1.0);
}
