#version 330

layout(location = 0) out vec4 fragColor;

uniform sampler2D font;

flat in int v_number;

void main() {

  vec2 st = gl_PointCoord;
  st.x /= 10.;
  st.x += (v_number - 1) / 10.0;

  vec4 color = texture(font, st);
  if (color.a < 1e-6) discard;
  fragColor = vec4(1, 0, 0, 1);
}