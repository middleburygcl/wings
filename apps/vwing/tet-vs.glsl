#version 410
layout (location = 0 ) in vec3 a_Position;

flat out int v_id;

void main() {
  gl_Position = vec4(a_Position, 1.0);
  v_id = gl_VertexID;
}
