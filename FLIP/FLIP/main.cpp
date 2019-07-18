#include <glad\glad.h>
#include <GLFW\glfw3.h>

#include <iostream>
#include "CShader.h"
#include "CGLFW.h"

int main()
{
	GLFW ourGlfw;
	
	//Init & create window
	ourGlfw.setWindowTitle("FLIP");
	ourGlfw.init();

	//Create Shader
	Shader ourShader("./Shader/VertexShader.glsl", "./Shader/FragmentShader.glsl");
	unsigned int VBO, VAO;
	
	//draw buffer
	while (!ourGlfw.getWindowShouldClose())
	{
		ourGlfw.processInput();
		ourShader.use();
		ourGlfw.swapBuffers();
		Sleep(100);
	}

	ourGlfw.terminate();

	return 0;
}