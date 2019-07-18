#ifndef SAHDER_H
#define SHADER_H

#include <glad/glad.h>

#include <string>
#include <fstream>
#include <sstream>
#include <iostream>

class Shader
{
public:
	unsigned int ID;

	Shader(const GLchar* vertexPath, const GLchar* fragmentPath);
	Shader(const GLchar* vertexPath, const GLchar* fragmentPath, const GLchar* geometryPath);

	void use();
	void setBool(const std::string &name, bool value)const;
	void setInt(const std::string &name, int value)const;
	void setFloat(const std::string &name, float value)const;
	void setVec3(const std::string &name, float v0, float v1, float v2)const;
	void setVec3Array(const std::string &name, int count, float* value)const;

private:
	void checkCompileError(unsigned int Shader, std::string type);
};
#endif //!SHADER_H