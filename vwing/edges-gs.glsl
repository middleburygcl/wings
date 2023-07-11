#version 330 core

layout (points) in;
uniform mat4 u_ModelViewProjectionMatrix;
uniform mat4 u_NormalMatrix;
uniform mat4 u_ModelViewMatrix;
uniform mat4 u_ViewportMatrix;

uniform int u_clip;
uniform vec3 u_clip_point;
uniform vec3 u_clip_normal;

uniform vec2 u_ViewportSize = vec2(800,600);
uniform vec2 u_aa_radius;

flat out int id;

//layout (line_strip , max_vertices = 3) out;
layout (triangle_strip, max_vertices = 6) out;

float t = 0.025;

uniform samplerBuffer coordinates;
uniform usamplerBuffer index;

flat in int v_id[];

void main() {

  uint i0 = texelFetch(index, 2 * v_id[0]    ).r;
  uint i1 = texelFetch(index, 2 * v_id[0] + 1).r;

  vec3 x0 = texelFetch(coordinates, int(i0)).xyz;
  vec3 x1 = texelFetch(coordinates, int(i1)).xyz;

  vec4 p0 = u_ModelViewProjectionMatrix * vec4(x0, 1.0);
  vec4 p1 = u_ModelViewProjectionMatrix * vec4(x1, 1.0);

#if 0
  gl_Position = p0;
  id = v_id[0];//gl_PrimitiveIDIn;
  EmitVertex();

  gl_Position = p1;
  id = v_id[0];//gl_PrimitiveIDIn;
  EmitVertex();

  EndPrimitive();
#elif 0

  float len = length(p1.xyz - p0.xyz);
  vec3 dir = normalize(p1.xyz - p0.xyz);
  vec3 n = vec3(dir[1], -dir[0], 0);

  float w = 0.008 * max(abs(p0.w), abs(p1.w));

  vec3 q0 = p0.xyz - 0.000 * len * dir - w * n;
  vec3 q1 = p1.xyz + 0.000 * len * dir - w * n;
  vec3 q2 = p1.xyz + 0.000 * len * dir + w * n;
  vec3 q3 = p0.xyz - 0.000 * len * dir + w * n;

  gl_Position = vec4(q0, p0.w);
  id = v_id[0];
  EmitVertex();

  gl_Position = vec4(q1, p1.w);
  id = v_id[0];
  EmitVertex();

  gl_Position = vec4(q2, p1.w);
  id = v_id[0];
  EmitVertex();

  EndPrimitive();

  gl_Position = vec4(q0, p0.w);
  id = v_id[0];
  EmitVertex();

  gl_Position = vec4(q2, p1.w);
  id = v_id[0];
  EmitVertex();

  gl_Position = vec4(q3, p0.w);
  id = v_id[0];
  EmitVertex();

  EndPrimitive();

#else
  float u_width        = u_ViewportSize[0];
  float u_height       = u_ViewportSize[1];
  float u_aspect_ratio = u_height / u_width;

  vec2 v_line_width;
  v_line_width[0] = 5;
  v_line_width[1] = 5;

  vec2 ndc_a = p0.xy / p0.w;
  vec2 ndc_b = p1.xy / p1.w;

  vec2 line_vector = ndc_b - ndc_a;
  vec2 viewport_line_vector = line_vector * u_ViewportSize;
  vec2 dir = normalize(vec2(line_vector.x, line_vector.y * u_aspect_ratio));
  
  float line_width_a     = max( 1.0, v_line_width[0] ) + u_aa_radius[0];
  float line_width_b     = max( 1.0, v_line_width[1] ) + u_aa_radius[0];
  float extension_length = u_aa_radius[1];
  float line_length      = length( viewport_line_vector ) + 2.0 * extension_length;
  
  vec2 normal    = vec2( -dir.y, dir.x );
  vec2 normal_a  = vec2( line_width_a/u_width, line_width_a/u_height ) * normal;
  vec2 normal_b  = vec2( line_width_b/u_width, line_width_b/u_height ) * normal;
  vec2 extension = vec2( extension_length / u_width, extension_length / u_height ) * dir;


  //g_col = vec4( v_col[0].rgb, v_col[0].a * min( v_line_width[0], 1.0f ) );
  //g_u = line_width_a;
  //g_v = line_length * 0.5;
  //g_line_width = line_width_a;
  //g_line_length = line_length * 0.5;
  gl_Position = vec4( (ndc_a + normal_a - extension) * p0.w, p0.zw );
  id = v_id[0];
  EmitVertex();
  
  //g_u = -line_width_a;
  //g_v = line_length * 0.5;
  //g_line_width = line_width_a;
  //g_line_length = line_length * 0.5;
  gl_Position = vec4( (ndc_a - normal_a - extension) * p0.w, p0.zw );
  id = v_id[0];
  EmitVertex();
  
  //g_col = vec4( v_col[0].rgb, v_col[0].a * min( v_line_width[0], 1.0f ) );
  //g_u = line_width_b;
  //g_v = -line_length * 0.5;
  //g_line_width = line_width_b;
  //g_line_length = line_length * 0.5;
  gl_Position = vec4( (ndc_b + normal_b + extension) * p1.w, p1.zw );
  id = v_id[0];
  EmitVertex();
  
  //g_u = -line_width_b;
  //g_v = -line_length * 0.5;
  //g_line_width = line_width_b;
  //g_line_length = line_length * 0.5;
  gl_Position = vec4( (ndc_b - normal_b + extension) * p1.w, p1.zw );
  id = v_id[0];
  EmitVertex();
  
  EndPrimitive();
#endif
}
