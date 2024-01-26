#version 330 core

layout (points) in;
layout (points, max_vertices = 8) out;

uniform mat4 u_ModelViewProjectionMatrix;

flat out int v_number;

uniform samplerBuffer points;
uniform usamplerBuffer visibility;

// TODO pass uniforms
float width = 800;
float height = 600;

flat in int v_id[];

int output_character(in float px, in vec4 q, in int n) {
  gl_Position = q;
  gl_PointSize = 20;
  v_number = int(mod(n, 10));
  EmitVertex();

  if (n < 10) return -1;
  return n / 10;
}

void main() {

  // retrieve point number and initial point coordinates
  int n = v_id[0];
  int visible = int(texelFetch(visibility, n).r);
  if (visible == 0) return;
  vec3 p = texelFetch(points, n).xyz;
  vec4 q = u_ModelViewProjectionMatrix * vec4(p * 1.00001, 1);

  float w = q.w; // depth
  float px = ((q.x + 1) * width / 2) / w; // pixel x-coordinate

  n = output_character(px, q, n);
  while (n >= 0) {
    px -= 22; // each character in the bitmap has a width of 22 pixels
    q.x = 2 * w * px / width - 1;
    n = output_character(px, q, n);
  }

  EndPrimitive();
}