#version 410 core
// LEOPOLDO ZUGASTI 260919951
in vec3 vertex;
in vec3 normal;
in float ut;
in float phi;

out vec3 camSpacePosition;
out vec3 camSpaceNormal;
out float utv;
out float phiv;

uniform mat4 P;
uniform mat4 MV;

uniform vec3 lightCamSpacePosition;
uniform vec3 lightColor;
uniform vec3 materialDiffuse;
uniform float materialShininess; 

void main() {	
	vec4 csp = MV * vec4( vertex, 1 );
	camSpacePosition = csp.xyz;
	utv = ut;
	phiv = phi;
	// assuming view is a rigid transform, just use the view matrix for the normal
	camSpaceNormal = normalize((MV * vec4( normal, 0 )).xyz);	
	gl_Position = P * csp;
}
