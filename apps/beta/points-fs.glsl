#version 330

layout(location = 0) out vec4 fragColor;

void main() {
  float alpha = 1 / (1 + pow(5 * length(gl_PointCoord - vec2(0.5, 0.5)), 4));
  if (alpha < 0.5) discard;
  fragColor = vec4(0, 0, 0, alpha);
}
