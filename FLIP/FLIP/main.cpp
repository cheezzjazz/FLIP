#include <glad\glad.h>
#include <GLFW\glfw3.h>

#include <iostream>
#include <vector>
#include "CShader.h"
#include "CGLFW.h"
#include "CSolver.h"

#define PARTICLE_COUNTS 10
#define GRID_SIZE 20

int main()
{
	GLFW ourGlfw;
	
	//Init & create window
	ourGlfw.setWindowTitle("FLIP 2D");
	ourGlfw.init();
	
	float* vertices = new float[PARTICLE_COUNTS * 2];	//2 means x, y
	for (int i = 0; i < PARTICLE_COUNTS; i++)
	{
		vertices[2 * i] = 0.8f - 0.2f*i;
		vertices[2 * i + 1] = 0.5f  - 0.2f*i;
	}

	//FLIP INIT
	Solver ourSolver(GRID_SIZE);
	

	//Create Shader
	Shader ourShader("./Shader/VertexShader.glsl", "./Shader/FragmentShader.glsl");
	unsigned int VBO, VAO;
	glGenBuffers(1, &VBO);
	glGenVertexArrays(1, &VAO);
	
	glBindVertexArray(VAO);
	glBindBuffer(GL_ARRAY_BUFFER, VBO);
	glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);

	glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 2 * sizeof(float), (void*)0);
	glEnableVertexAttribArray(0);
	
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

		glBindBuffer(GL_ARRAY_BUFFER, VBO);
		float * ptr = (float*) glMapBuffer(GL_ARRAY_BUFFER, GL_READ_WRITE);
		if (ptr)
		{
			ourSolver.stepFLIP();
			glUnmapBuffer(GL_ARRAY_BUFFER);
		}

		glPointSize(5.0f);
		//glDrawArrays(GL_POINTS, 0, PARTICLE_COUNTS);
	
		ourGlfw.swapBuffers();
		
		Sleep(100);
		--nFrame;
	}

	glDeleteBuffers(1, &VBO);
	glDeleteVertexArrays(1, &VAO);

	ourGlfw.terminate();

	return 0;
}