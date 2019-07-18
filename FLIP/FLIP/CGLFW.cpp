#include "CGLFW.h"

static void errorCallback(int errorCode, const char* errorDescription)
{
	std::cout << "ERROR::GLFW " << errorDescription << std::endl;
}

void framebuffer_size_callback(GLFWwindow* window, int width, int height)
{
	glViewport(0, 0, width, height);
}

void GLFW::setWindowTitle(char * title)
{
	Window_Title = title;
}

void GLFW::loadGLADLoader()
{
	if (!gladLoadGLLoader((GLADloadproc)(glfwGetProcAddress)))
	{
		std::cout << "Failed to initialize GLAD" << std::endl;
		glfwTerminate();
		return;
	}
}

int GLFW::init()
{
	glfwSetErrorCallback(errorCallback);
	if (!glfwInit())
	{
		std::cout << "ERROR::GLFW:: Failed to initialize GLFW" << std::endl;
		return -1;
	}
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

	window = glfwCreateWindow(scr_width, scr_height, Window_Title, NULL, NULL);
	if (window == NULL)
	{
		std::cout << "Failed to create GLFW window " << std::endl;
		glfwTerminate();

		return -1;
	}

	glfwMakeContextCurrent(window);
	glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);

	loadGLADLoader();
	glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT);
	return 0;
}

void GLFW::swapBuffers()
{
	if (!checkWindowcreated())
		return;

	glfwSwapBuffers(window);
	glfwPollEvents();
}

bool GLFW::getKeyPressed(const int & key)
{
	if (!checkWindowcreated())
		return 0;

	return (glfwGetKey(window, key) == GLFW_PRESS);
}

int GLFW::checkWindowcreated()
{
	if (window == NULL)
	{
		std::cout << "ERROR::GLFW:: Window was not created yet!" << std::endl;
		return 0;
	}
	return 1;
}

bool GLFW::getWindowShouldClose()
{
	if (!checkWindowcreated())
		return 0;

	return (glfwWindowShouldClose(window) != 0);
}

void GLFW::terminate()
{
	glfwTerminate();
}

void GLFW::processInput()
{
	if (getKeyPressed(GLFW_KEY_ESCAPE))
		glfwSetWindowShouldClose(window, true);
	else if (getKeyPressed(GLFW_KEY_SPACE))
		isStart = true;
}

