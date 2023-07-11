#version 330
layout(location = 0) in vec3 a_Position;

uniform mat4 u_ModelViewProjectionMatrix;

void main() {
#if WITH_GS
  gl_Position = vec4(a_Position, 1.0);
#else
  gl_Position = u_ModelViewProjectionMatrix * vec4(a_Position, 1.0);
  gl_PointSize = 20.0;
#endif
}
