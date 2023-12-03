#version 410 core
// LEOPOLDO ZUGASTI 260919951
in vec3 fragColor;
out vec4 out_fragColor;

void main()
{
    out_fragColor = vec4(fragColor, 1.0);
}
