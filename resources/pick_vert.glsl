#version 410 core
// LEOPOLDO ZUGASTI 260919951
in vec3 vertex;
in vec3 color;

out vec3 colorv;

uniform mat4 P;
uniform mat4 MV;

void main() {
	colorv = color;	
	gl_Position = P * MV * vec4( vertex, 1 );
}