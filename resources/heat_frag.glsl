#version 410 core

// LEOPOLDO ZUGASTI 260919951

in vec3 camSpacePosition;
in vec3 camSpaceNormal;
in float utv;
in float phiv;

out vec4 out_fragColor;

uniform vec3 lightCamSpacePosition;
uniform vec3 lightColor;
uniform vec3 materialDiffuse;
uniform float materialShininess; 

void main(void) {
	
	vec3 v = normalize(-camSpacePosition);
	vec3 n = normalize(camSpaceNormal);
	vec3 l = normalize(lightCamSpacePosition - camSpacePosition);

	// TODO: 4, 7, 11 Implement your GLSL per fragement lighting, heat colouring, and distance stripes here!

	vec3 h = normalize(l + v);
	float specular = max(dot( n, h ), 0);
	float diffuse = max(dot(n, l), 0);

	specular = pow(specular, materialShininess);
	vec3 hotColor = vec3(1.0, 0.0, 0.0);
	vec3 neutralColor = vec3(0.5, 0.5, 0.5); // This is purple (between blue and red)
	vec3 newColor = vec3(0.0, 0.0, 0.0);
	float newUtv = (utv + 0.5) / 1.5;
	float stripeDetection = mod(phiv, 1.0/150.0);//cos(2*3.1415926* 14*(phiv));//mod(smoothstep(-1.0,1.0,phiv)*(1.0-smoothstep(1.0,-1.0,phiv)), 1/150.0);
	// if we are on a distance stripe
	if (phiv != 0 && stripeDetection < 0.001 && stripeDetection > -0.001)
	{
		newColor = vec3(1.0, 1.0, 1.0);
	}
	// if we are not on a stripe
	else 
	{
		if (newUtv <= 0.5)
		{
			newUtv = clamp(newUtv * 2.0, 0, 1);
			newColor = (1 - newUtv) * materialDiffuse + (newUtv) * neutralColor;
		}
		else
		{
			newUtv = clamp((newUtv - 0.5)* 2.0, 0, 1);
			newColor = (1 - newUtv) * neutralColor + (newUtv) * hotColor;
		}
	}
	
	vec3 reflectedLight = newColor * lightColor * specular;
	vec3 ambientLight = newColor * vec3(0.1, 0.1, 0.1);
	vec3 scatteredLight = (newColor/2.0) * lightColor * diffuse;

	out_fragColor = vec4( reflectedLight + ambientLight + scatteredLight, clamp(n.x + utv + phiv, 0, 0) + 1.0 );
}
