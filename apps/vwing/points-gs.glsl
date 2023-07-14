#version 330 core

layout (points) in;
uniform mat4 u_ModelViewProjectionMatrix;

layout (triangle_strip , max_vertices = 60) out;

uniform float u_length;
const float t = (1.0 + sqrt(5.0)) / 2.0;
const vec3 pts[12] = vec3[12](vec3(t, 1, 0), vec3(-t, 1, 0), vec3(t, -1, 0), vec3(-t, -1, 0),
                      vec3(1, 0, t), vec3(1, 0, -t), vec3(-1, 0, t), vec3(-1, 0, -t),
                      vec3(0, t, 1), vec3(0, -t, 1), vec3(0, t, -1), vec3(0, -t, -1));

const ivec3 tris[20] = ivec3[20](
      ivec3(0, 8, 4),   // 0
      ivec3(0, 5, 10),  // 1
      ivec3(2, 4, 9),   // 2
      ivec3(2, 11, 5),  // 3
      ivec3(1, 6, 8),   // 4
      ivec3(1, 10, 7),  // 5
      ivec3(3, 9, 6),   // 6
      ivec3(3, 7, 11),  // 7
      ivec3(0, 10, 8),  // 8
      ivec3(1, 8, 10),  // 9
      ivec3(2, 9, 11),  // 10
      ivec3(3, 11, 9),  // 11
      ivec3(4, 2, 0),   // 12
      ivec3(5, 0, 2),   // 13
      ivec3(6, 1, 3),   // 14
      ivec3(7, 3, 1),   // 15
      ivec3(8, 6, 4),   // 16
      ivec3(9, 4, 6),   // 17
      ivec3(10, 5, 7),  // 18
      ivec3(11, 7, 5)   // 19
);

void make_triangle(vec3 x0 , vec3 x1 , vec3 x2) {

  vec4 p0 = u_ModelViewProjectionMatrix * vec4(x0, 1.0);
  vec4 p1 = u_ModelViewProjectionMatrix * vec4(x1, 1.0);
  vec4 p2 = u_ModelViewProjectionMatrix * vec4(x2, 1.0);

  gl_Position = p0;
  EmitVertex();

  gl_Position = p1;
  EmitVertex();

  gl_Position = p2;
  gl_PrimitiveID = gl_PrimitiveIDIn;
  EmitVertex();

  EndPrimitive();
}

void main() {
  float r = 0.01 * u_length;

  vec3 c = gl_in[0].gl_Position.xyz;
  for (int i = 0; i < 20; i++) {
    int p0 = tris[i].x;
    int p1 = tris[i].y;
    int p2 = tris[i].z;
    make_triangle(c + r * pts[p0], c + r * pts[p1], c + r * pts[p2]);
  }
}
