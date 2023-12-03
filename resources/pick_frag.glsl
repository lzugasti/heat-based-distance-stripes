#version 410 core
// LEOPOLDO ZUGASTI 260919951
in vec3 colorv;
out vec4 out_fragColor;

void main(void) {
	out_fragColor = vec4( colorv / 255.0f, 1 );
}