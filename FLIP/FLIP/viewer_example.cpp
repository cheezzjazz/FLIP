#include <glad\glad.h>
#include <GLFW\glfw3.h>

#include <iostream>
#include "CShader.h"
#include "CGLFW.h"
#include "CSolver.h"

int main()
{
	GLFW ourGlfw;

	//Init & create window
	ourGlfw.setWindowTitle("FLIP 2D");
	ourGlfw.init();

	float vertices[] =
	{
		-0.5f, 0.5f, 1.0f, 1.0f, 1.0f,// top mid 0
		0.0f, 0.5f,  1.0f, 1.0f, 1.0f,// bottom right 1
		0.0f, -0.5f, 1.0f, 1.0f, 1.0f,// // bottom left 2
		-0.5f, -0.5f, 1.0f, 1.0f, 1.0f,

		-0.25f, 0.3f, 0.0f, 0.0f, 1.0f,// top mid 0
		0.5f, 0.3f, 0.0f, 0.0f, 1.0f,// bottom right 1
		0.5f, -0.8f,0.0f, 0.0f, 1.0f,// // bottom left 2
		-0.25f, -0.8f, 0.0f, 0.0f, 1.0f
	};

	unsigned int indices[] = {	//note that we start from 0!!
		0, 1, 2,	// triangle
		2, 3, 0,
		4, 5, 6,
		6, 7, 4
	};
	//Create Shader
	Shader ourShader("./Shader/VertexShader_color.glsl", "./Shader/FragmentShader.glsl");
	unsigned int VBO, VAO, EBO;
	glGenBuffers(1, &VBO);
	glGenVertexArrays(1, &VAO);

	glGenBuffers(1, &EBO);

	glBindVertexArray(VAO);
	glBindBuffer(GL_ARRAY_BUFFER, VBO);
	glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);

	//glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 2 * sizeof(float), (void*)0);
	glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 5 * sizeof(float), (void*)0);
	glEnableVertexAttribArray(0);
	//color
	glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 5 * sizeof(float), (void*)(3 * sizeof(float)));
	glEnableVertexAttribArray(1);

	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(indices), indices, GL_STATIC_DRAW);

	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindVertexArray(0);

	glBindVertexArray(VAO);
	//draw buffer

	int nFrame = 500;
	while (!ourGlfw.getWindowShouldClose() || nFrame == 0)
	{
		std::cout << "Frame [" << nFrame << "] \n";

		ourGlfw.processInput();
		ourShader.use();

		//glPointSize(2.0f);
		//glDrawArrays(GL_POINTS, 0, 8);
		glDrawElements(GL_TRIANGLES, 12, GL_UNSIGNED_INT, 0);
		ourGlfw.swapBuffers();

		Sleep(100);
		--nFrame;
	}

	glDeleteBuffers(1, &VBO);
	glDeleteVertexArrays(1, &VAO);

	ourGlfw.terminate();

	return 0;
}