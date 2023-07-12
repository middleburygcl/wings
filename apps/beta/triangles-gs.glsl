#version 330 core

layout (points) in;
uniform mat4 u_ModelViewProjectionMatrix;
uniform mat4 u_NormalMatrix;
uniform mat4 u_ModelViewMatrix;
uniform mat4 u_ViewportMatrix;

uniform int u_clip;
uniform vec3 u_clip_point;
uniform vec3 u_clip_normal;

uniform vec2 u_ViewportSize;// = vec2(800,600);

noperspective out vec3 v_Position;
noperspective out vec3 altitude;
out vec3 v_Normal;
noperspective out vec2 v_Barycentric;
flat out int id;

layout (triangle_strip , max_vertices = 3) out;

uniform samplerBuffer coordinates;
uniform usamplerBuffer index;

flat in int v_id[];

void main() {

  uint i0 = texelFetch(index, 3 * v_id[0]    ).r;
  uint i1 = texelFetch(index, 3 * v_id[0] + 1).r;
  uint i2 = texelFetch(index, 3 * v_id[0] + 2).r;

  vec3 x0 = texelFetch(coordinates, int(i0)).xyz;
  vec3 x1 = texelFetch(coordinates, int(i1)).xyz;
  vec3 x2 = texelFetch(coordinates, int(i2)).xyz;

  vec4 u = vec4(x1 - x0, 0.0);
  vec4 v = vec4(x2 - x0, 0.0);
  vec3 n = normalize(mat3(u_NormalMatrix) * cross(u.xyz, v.xyz));

  if (u_clip > 0) {
    float result =
      sign(dot(u_clip_normal, x0 - u_clip_point)) +
      sign(dot(u_clip_normal, x1 - u_clip_point)) +
      sign(dot(u_clip_normal, x2 - u_clip_point));
    if (result < 3) return;
  }

  vec4 p0 = u_ModelViewProjectionMatrix * vec4(x0, 1.0);
  vec4 p1 = u_ModelViewProjectionMatrix * vec4(x1, 1.0);
  vec4 p2 = u_ModelViewProjectionMatrix * vec4(x2, 1.0);

  vec2 q0 = p0.xy / p0.w;
  vec2 q1 = p1.xy / p1.w;
  vec2 q2 = p2.xy / p2.w;

  vec2 v1 = u_ViewportSize * (q1 - q0);
  vec2 v2 = u_ViewportSize * (q2 - q0);
  vec2 v3 = u_ViewportSize * (q2 - q1);

  float area = abs(v1.x * v2.y - v1.y * v2.x); // twice the area
  float h1 = area / length(v3);
  float h2 = area / length(v2);
  float h3 = area / length(v1);

  gl_Position = p0;
  v_Position  = (u_ModelViewMatrix * vec4(x0, 1.0)).xyz;
  altitude    = vec3(h1, 0, 0);
  id = v_id[0];//gl_PrimitiveIDIn;
  v_Normal = n;
  v_Barycentric = vec2(0,0);
  EmitVertex();

  gl_Position = p1;
  v_Position  = (u_ModelViewMatrix * vec4(x1, 1.0)).xyz;
  altitude    = vec3(0, h2, 0);
  id = v_id[0];//gl_PrimitiveIDIn;
  v_Normal = n;
  v_Barycentric = vec2(1,0);
  EmitVertex();

  gl_Position = p2;
  v_Position  = (u_ModelViewMatrix * vec4(x2, 1.0)).xyz;
  altitude    = vec3(0, 0, h3);
  id = v_id[0];//gl_PrimitiveIDIn;
  v_Normal = n;
  v_Barycentric = vec2(0,1);

  gl_PrimitiveID = gl_PrimitiveIDIn;

  EmitVertex();

  EndPrimitive();
}
