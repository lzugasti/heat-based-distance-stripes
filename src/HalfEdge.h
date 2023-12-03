// LEOPOLDO ZUGASTI 260919951

#pragma once
#ifndef HalfEdge_H
#define HalfEdge_H
#include <memory>

#include "GLHeaders.h"

#include "Vertex.h"
#include "Program.h"
#include "MatrixStack.h"

using namespace std;

class Face;

/**
 * Class containing the half edge pointers, and a method for drawing
 * the half edge for debugging and evaluation.
 */
class HalfEdge
{
public:
	/* pointer to twin of this half edge */
	HalfEdge *twin;

	/* pointer to next half edge */
	HalfEdge *next;

	/* pointer to head of half edge */
	Vertex *head;

	/* pointer to face at the left side */
	Face *leftFace;

	/**
	 * edge opposite head of this half edge cross the normal
	 * (only makes sense for triangles)
	 */
	glm::vec3 ecn;

	/**
	 * this edge vector
	 */
	glm::vec3 e;

	/*	@return pointer to the previous half edge (could just be stored to save time) */
	HalfEdge *prev();

	/**
	 * Displays the half edge as a half arrow pointing to the head vertex.
	 */
	void display(shared_ptr<Program> prog, shared_ptr<MatrixStack> P, shared_ptr<MatrixStack> MV, GLuint vao, GLuint vbo);
};
#endif
