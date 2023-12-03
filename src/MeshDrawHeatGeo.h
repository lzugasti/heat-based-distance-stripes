// LEOPOLDO ZUGASTI 260919951

#pragma once
#ifndef MeshDrawHeatGeo_H
#define MeshDrawHeatGeo_H

#include "GLHeaders.h"

#include "GLSL.h"
#include "HEDS.h"
#include "Program.h"
#include "MatrixStack.h"
/**
 * Class for organizing GLSL programs and mesh data for picking triangles and
 * likewise vertices based on barycentric coordinates
 * @author kry
 */
class MeshDrawHeatGeo
{
public:
	/* height of the gradient from the surface */
	float gradHeight = 0.1f; // [1e-3, 1e2]
	/* scale of the gradient */
	float gradScale = 0.5; // [1e-3, 1e3]

	/* scale distance solution values */
	float scale = 0.5f; // [1e-4, 1e4]

	// TODO: feel free to change the value above for better visualization!

	shared_ptr<HEDS> heds = NULL;

	/** Vertex and normal data */
	GLuint pvdVBOID = 0;
	/** Dynamic vertex data */
	GLuint dvdVBOID = 0;
	/** Index buffer for drawing */
	GLuint indexVBOID = 0;
	/** Gradients */
	GLuint gradientBufferID = 0;

	int pvdCount = 0;

	int pvdStrideBytes = 0;
	int dvdStrideBytes = 0;

	int gradientBufferSize = 0;

	GLfloat *pvdBuffer = NULL;
	GLfloat *dvdBuffer = NULL;
	GLuint *indexBuffer = NULL;
	GLfloat *gradientBuffer = NULL;

	shared_ptr<Program> aProgram = NULL;

	// vao
	GLuint vao = 0;
	// gradient vao
	GLuint gradientVao = 0;

	// VBO ID for camSpacePosition
	int camSpacePositionID = 0;
	// VBO ID for camSpaceNormal
	int camSpaceNormalID = 0;
	// VBO ID for utv
	int utvID = 0;
	// VBO ID for phiv
	int phivID = 0;

	/**
	 * Initializes the GLSL shader (the mesh must also be set before picking)
	 *
	 * @param pProgram default drawing program
	 */
	void init(shared_ptr<Program> pProgram);

	/**
	 * Load a new mesh
	 *
	 * @param heds new mesh
	 */
	void setMesh(shared_ptr<HEDS> heds);

	/**
	 * Fill dynamic VBO infos
	 */
	void fillDynamicVBO();

	/**
	 * fill static VBO infos
	 */
	void fillStaticVBOs();

	/**
	 * Draw VBOs using the default shader program
	 *
	 * @param P projection
	 * @param MV modelview
	 */
	void drawVBOs(shared_ptr<MatrixStack> P, shared_ptr<MatrixStack> MV);

	/**
	 * Draw gradients using the provided shader program
	 *
	 * @param P projection
	 * @param MV modelview
	 * @param lineProg program for line drawing
	 */
	void drawGrad(shared_ptr<MatrixStack> P, shared_ptr<MatrixStack> MV, shared_ptr<Program> lineProg);
};
#endif