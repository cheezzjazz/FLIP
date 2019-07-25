#include <glad\glad.h>
#include <GLFW\glfw3.h>

#include <iostream>
#include <vector>
#include "CShader.h"
#include "CGLFW.h"
#include "CSolver.h"

#define GRID_SIZE 128
void updateVertices(float* mappedbuffer, float* srcvertices, int count)
{
	if (!mappedbuffer || !srcvertices)
	{
		std::cout << "ERROR::UPDATEVERTICES:: DATA Pointer is NULL" << std::endl;
		return;
	}
	float diff = 0.01f;
	//diff += 0.1f;
	std::cout << diff << std::endl;
	for (int i = 0; i < count; i++)
	{
		//float x, y;
		//x = *srcvertices; srcvertices++;
		//y = *srcvertices; srcvertices++;
		//std::cout << "x=" << x << "y=" << y;
		////float diff = 0.008f;
		//mappedbuffer[2*i] = x + diff;
		//srcvertices[2 * i] = mappedbuffer[2 * i];
		////mappedbuffer++;
		//mappedbuffer[2*i+1] = y + diff;
		//srcvertices[2 * i] = mappedbuffer[2 * i];
		//std::cout << "map[" << 2 * i << "]=" << mappedbuffer[2 * i] << "map[" << 2 * i + 1 << "]=" << mappedbuffer[2 * i + 1] << std::endl;
		////mappedbuffer++;
		////update vertex coords
		*srcvertices += diff;
		*mappedbuffer = *srcvertices; std::cout << "x:" << *mappedbuffer;
		++srcvertices; ++mappedbuffer;
		*srcvertices += diff;
		*mappedbuffer = *srcvertices + diff; std::cout << "y:" << *mappedbuffer<<std::endl;
		++srcvertices; ++mappedbuffer;
		//*srcvertices += diff;
		//*mappedbuffer = *srcvertices; std::cout << "z:" << *mappedbuffer << std::endl;
		//++srcvertices; ++mappedbuffer;
	}


}

int main()
{
	GLFW ourGlfw;
	
	//Init & create window
	ourGlfw.setWindowTitle("FLIP 2D");
	ourGlfw.init();
	
	//FLIP INIT
	Solver ourSolver(GRID_SIZE);
	int particlecnt = 0;
	particlecnt = ourSolver.m_particles->num;
	std::cout << "frame : particle count : " << particlecnt << std::endl;
	float* vertices = new float[particlecnt * 2];	//2 means x, y
	for (int i = 0; i < particlecnt; i++)
	{
		vertices[2 * i] = 0.0f + 0 ;
		vertices[2 * i + 1] = 0.0f;
	}


	//Create Shader
	Shader ourShader("./Shader/VertexShader.glsl", "./Shader/FragmentShader.glsl");
	unsigned int VBO, VAO;
	glGenBuffers(1, &VBO);
	glGenVertexArrays(1, &VAO);
	
	glBindVertexArray(VAO);
	glBindBuffer(GL_ARRAY_BUFFER, VBO);
	glBufferData(GL_ARRAY_BUFFER, 2*particlecnt*sizeof(float), vertices, GL_STATIC_DRAW);

	glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 2 * sizeof(float), (void*)0);
	glEnableVertexAttribArray(0);
	
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindVertexArray(0);

	glBindVertexArray(VAO);
	
	//draw buffer
	int nFrame = 100;
	while (!ourGlfw.getWindowShouldClose())// && nFrame != 0)
	{
		glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
		glClear(GL_COLOR_BUFFER_BIT);
		//std::cout << "Frame [" << nFrame << "] \n";

		ourGlfw.processInput();
		ourShader.use();

		glBindBuffer(GL_ARRAY_BUFFER, VBO);
		float * ptr = (float*) glMapBuffer(GL_ARRAY_BUFFER, GL_READ_WRITE);
		if (ptr)
		{
			//ourSolver.stepFLIP();
			//ourSolver.update(ptr, vertices);
			ourSolver.update(ptr, vertices, ourGlfw.scr_width, ourGlfw.scr_height);
			//updateVertices(ptr, vertices, particlecnt);
			glUnmapBuffer(GL_ARRAY_BUFFER);
		}

		glPointSize(2.0f);
		glDrawArrays(GL_POINTS, 0, particlecnt);
	
		ourGlfw.swapBuffers();
		
		Sleep(100);
		--nFrame;
	}

	glDeleteBuffers(1, &VBO);
	glDeleteVertexArrays(1, &VAO);

	ourGlfw.terminate();

	return 0;
}