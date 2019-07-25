#version 330 core
layout (location = 0) in vec3 aPos;

out vec3 ourColor;

void main()
{
	gl_Position = vec4(aPos, 1.0f);
	ourColor = vec3(0.0f/255.0f, 128.0f/255.0f, 255.0f/255.0f);
}