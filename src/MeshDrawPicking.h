// LEOPOLDO ZUGASTI 260919951

#pragma once
#ifndef MeshDrawPicking_H
#define MeshDrawPicking_H

#include <memory>

#include "GLHeaders.h"

#include "GLSL.h"
#include "HEDS.h"
#include "Program.h"
/**
 * Class for organizing GLSL programs and mesh data for picking triangles and
 * likewise vertices based on barycentric coordinates
 * @author kry
 */
class MeshDrawPicking
{
public:
	/* Half edge data structure */
	shared_ptr<HEDS> heds = NULL;
	/* number of faces */
	int nFaces = 0;

	/* local vao */
	GLuint pickVao = 0;

	/* picked triangle VBO ID */
	GLuint pickedVBOID = 0;
	/* picked triangle buffer */
	vector<GLfloat> pickedTriangleVerticesBuffer;

	/* picking drawing vertex VBO ID */
	GLuint pickingVertexBOID = 0;
	/* picking drawing face ID VBO ID */
	GLuint pickingFaceIDBOID = 0;
	/* picking drawing triangle indexing VBO ID */
	GLuint pickingTriangleIndexBOID = 0;
	/* picking drawing vertex buffer */
	vector<GLfloat> pickingVertexB;
	/* picking drawing face ID buffer */
	vector<GLubyte> pickingFaceIDBuffer;
	/* picking drawing triangle indexing buffer */
	vector<GLuint> pickingTriangleIndexBuffer;
	/* bytes of vertex buffer */
	int pickingVertexBOStrideBytes = 0;
	/* bytes of face buffer */
	int pickingFaceIDBOStrideBytes = 0;

	/* default drawing program */
	shared_ptr<Program> aProgram = NULL;

	/* vertex attribute ID */
	int vertexID = 0;
	/* color attribute ID */
	int colorID = 0;

	/**
	 * Initializes the GLSL shader (the mesh must also be set before picking)
	 * @param gl
	 */
	void init(shared_ptr<Program> pProgram);

	/**
	 * Sets the mesh, allocates and sets VBOs for picking
	 * @param gl
	 * @param heds
	 */
	void setMesh(shared_ptr<HEDS> heds);

	/**
	 * Sets up the picking buffers by creating vertex data for each triangle with
	 * independent colours set for each vertex.
	 */
	void fillPickingBuffers();

	/**
	 * Draws for picking with a unique colour for every triangle
	 *
	 * @param p clicked point
	 * @param BCC picked point color
	 * @param window window object
	 * @param P projection
	 * @param MV modelview
	 * @param triProg picked triangle drawing program
	 * @return int
	 */
	int pickTriangles(glm::vec2 *p, glm::vec3 *BCC, GLFWwindow *window, shared_ptr<MatrixStack> P, shared_ptr<MatrixStack> MV, shared_ptr<Program> triProg);
};
#endif