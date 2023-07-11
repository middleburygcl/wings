#version 330 core

layout (points) in;

uniform mat4 u_ModelViewProjectionMatrix;
uniform mat4 u_NormalMatrix;
uniform mat4 u_ModelViewMatrix;
uniform mat4 u_ViewportMatrix;

uniform float u_alpha;

uniform vec2 u_ViewportSize = vec2(800, 600);

noperspective out vec3 v_Position;
noperspective out vec3 altitude;

uniform samplerBuffer coordinates;
uniform usamplerBuffer index;
uniform usamplerBuffer first;
uniform usamplerBuffer count;

layout (triangle_strip , max_vertices = 50) out;

flat in int[] v_id;
flat out int id;
out vec3 v_Normal;

#define LARGE_DISTANCE 1000000 // anything larger than max(width, height) would be fine

uniform int u_clip;
uniform vec3 u_clip_point;
uniform vec3 u_clip_normal;

void main() {

  uint n = texelFetch(count, v_id[0]).r;
  uint f = texelFetch(first, v_id[0]).r;

  // the following can be used for convex polygons, but the second version is better for general polygons
  #if 0
  uint i0 = texelFetch(index, int(f)).r;
  uint i1 = texelFetch(index, int(f + 1)).r;

  vec3 X0 = texelFetch(coordinates, int(i0) ).xyz;
  vec3 X1 = texelFetch(coordinates, int(i1) ).xyz;

  vec4 x0 = vec4(X0, 1.0);
  vec4 x1 = vec4(X1, 1.0);

  // to keep things simple, we will only clip if the first two points are on the wrong side of the clipping plane
  if (u_clip > 0) {
    float result = sign(dot(u_clip_normal, X0 - u_clip_point)) +
                   sign(dot(u_clip_normal, X1 - u_clip_point));
    //if (result < 2) return;
  }

  vec4 p0 = u_ModelViewProjectionMatrix * x0;
  vec4 p1 = u_ModelViewProjectionMatrix * x1;

  for (int i = 2; i < n; i++) {

    uint i2 = texelFetch(index, int(f + i)).r;
    vec4 x2 = vec4(texelFetch(coordinates, int(i2)).xyz, 1.0);
    vec4 p2 = u_ModelViewProjectionMatrix * x2;

    vec4 u = x1 - x0;
    vec4 v = x2 - x0;
    vec3 normal = normalize(mat3(u_NormalMatrix) * cross(u.xyz,v.xyz));
    //if (floor(u_alpha) * normal[2] < 0) return; // normal is facing away, this will be false if there is any transparency active

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

    // these area interior, so set the distance far away
    if (i == 2 && n != 3) {
      h2 = LARGE_DISTANCE;
    }
    else if (i == n-1 && n != 3) {
      h3 = LARGE_DISTANCE;
    }
    else {
      h2 = LARGE_DISTANCE;
      h3 = LARGE_DISTANCE;
    }

    gl_Position = p0;
    v_Position  = x0.xyz;
    altitude    = vec3(h1, 0, 0);
    id = v_id[0];
    v_Normal = normal;
    EmitVertex();

    gl_Position = p1;
    v_Position  = x1.xyz;
    altitude    = vec3(0, h2, 0);
    id = v_id[0];
    v_Normal = normal;
    EmitVertex();

    gl_Position = p2;
    v_Position  = x2.xyz;
    altitude    = vec3(0, 0, h3);
    id = v_id[0];
    v_Normal = normal;

    gl_PrimitiveID = gl_PrimitiveIDIn;
    EmitVertex();

    EndPrimitive();

    p1 = p2;
    x1 = x2;
  }
  #else

  // compute the centroid
  vec4 x0 = vec4(0,0,0,0);
  uint on_side = uint(0);
  for (uint i = uint(0); i < n; i++) {
    uint j = texelFetch(index, int(f + i)).r;
    vec3 X = texelFetch(coordinates, int(j) ).xyz;
    on_side += uint(sign(dot(u_clip_normal, X - u_clip_point)));
    x0 = x0 + vec4(X,0.0);
  }
  x0 = x0 / float(n);
  //x0 = normalize(x0);
  x0.w = 1.0;

  vec4 p0 = u_ModelViewProjectionMatrix * x0;

  if (u_clip > 0) {
    if (on_side != n) return;
  }

  for (uint i = uint(0); i < n; i++) {
    uint j = (i == uint(n - uint(1))) ? uint(0) : i + uint(1);

    uint i1 = texelFetch(index, int(f + i)).r;
    vec4 x1 = vec4(texelFetch(coordinates, int(i1)).xyz, 1.0);
    vec4 p1 = u_ModelViewProjectionMatrix * x1;

    uint i2 = texelFetch(index, int(f + j)).r;
    vec4 x2 = vec4(texelFetch(coordinates, int(i2)).xyz, 1.0);
    vec4 p2 = u_ModelViewProjectionMatrix * x2;

    vec4 u = x1 - x0;
    vec4 v = x2 - x0;
    vec3 normal = normalize(mat3(u_NormalMatrix) * cross(u.xyz,v.xyz));
    //if (floor(u_alpha) * normal[2] < 0) return; // normal is facing away, this will be false if there is any transparency active

    vec2 q0 = p0.xy / p0.w;
    vec2 q1 = p1.xy / p1.w;
    vec2 q2 = p2.xy / p2.w;

    vec2 v1 = u_ViewportSize * (q1 - q0);
    vec2 v2 = u_ViewportSize * (q2 - q0);
    vec2 v3 = u_ViewportSize * (q2 - q1);

    float area = abs(v1.x * v2.y - v1.y * v2.x); // twice the area

    float h1 = area / length(v3);
    float h2 = LARGE_DISTANCE;
    float h3 = LARGE_DISTANCE;

    gl_Position = p0;
    v_Position  = x0.xyz;
    altitude    = vec3(h1, 0, 0);
    id = v_id[0];
    v_Normal = normal;
    EmitVertex();

    gl_Position = p1;
    v_Position  = x1.xyz;
    altitude    = vec3(0, h2, 0);
    id = v_id[0];
    v_Normal = normal;
    EmitVertex();

    gl_Position = p2;
    v_Position  = x2.xyz;
    altitude    = vec3(0, 0, h3);
    id = v_id[0];
    v_Normal = normal;

    gl_PrimitiveID = gl_PrimitiveIDIn;
    EmitVertex();

    EndPrimitive();
  }
  #endif

}
