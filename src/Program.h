// LEOPOLDO ZUGASTI 260919951

#pragma once
#ifndef Program_H
#define Program_H

#include <map>
#include <string>

#include "GLHeaders.h"

/**
 * An OpenGL Program (vertex and fragment shaders)
 */
class Program
{
public:
	Program();
	virtual ~Program();

	void setVerbose(bool v) { this->verbose = v; }

	void setShaderNames(const std::string &v, const std::string &f);
	virtual bool init();
	virtual void bind();
	virtual void unbind();

	void addAttribute(const std::string &name);
	void addUniform(const std::string &name);
	GLint getAttribute(const std::string &name) const;
	GLint getUniform(const std::string &name) const;

protected:
	std::string vShaderName;
	std::string fShaderName;

private:
	GLuint pid;
	std::map<std::string, GLint> attributes;
	std::map<std::string, GLint> uniforms;
	bool verbose = false;
};

#endif
