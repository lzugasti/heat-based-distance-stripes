// LEOPOLDO ZUGASTI 260919951

#pragma once
#ifndef HEDS_H
#define HEDS_H

#include <vector>
#include <map>
#include <memory>
#include <math.h>
#include "Face.h"
#include "PolygonSoup.h"

using namespace std;

/**
 * Half edge data structure.
 * Maintains a list of faces (i.e., one half edge of each) to allow
 * for easy display of geometry.
 */
class HEDS
{
public:
	~HEDS() {
		faces->clear();
		faces = NULL;
		vertices->clear();
		vertices = NULL;
		halfEdges->clear();
		halfEdges = NULL;
	}

	/**
	 * List of faces
	 */
	shared_ptr<vector<Face *>> faces = make_shared<vector<Face *>>();

	/**
	 * List of vertices
	 */
	shared_ptr<vector<Vertex *>> vertices = NULL;

	/** min phi value */
	double minphi = 0;

	/** max phi value */
	double maxphi = 0;

	/**
	 * Builds a half edge data structure from the polygon soup
	 * @param soup
	 */
	HEDS(shared_ptr<PolygonSoup> soup);

	/**
	 * Convenience map for keeping track of half edges you make or need
	 */
	shared_ptr<map<string, HalfEdge *>> halfEdges = make_shared<map<string, HalfEdge *>>();

	/**
	 * Helper function for creating a half edge, and pairing it up with its twin if
	 * the twin half edge has already been created.
	 *
	 * @param soup a polygon soup
	 * @param i index of the first vertex
	 * @param j index of the second vertex
	 * @return HalfEdge* pointer to the new HalfEdge object
	 */
	HalfEdge *createHalfEdge(shared_ptr<PolygonSoup> soup, int i, int j);

	/// The rest of this is about heat flow

	/**
	 * Initialize / reset heat flow
	 */
	void initHeatFlow();

	/**
	 * Perform a specified number of projected Gauss-Seidel steps of the heat
	 * diffusion equation. The current ut values stored in the vertices will be
	 * refined by this call.
	 *
	 * @param GSSteps number of steps to take
	 * @param t       solution time
	 */
	void solveHeatFlowStep(int GSSteps, double t);

	/**
	 * Precompute quantities needed to compute grad U, and div X
	 */
	void precomputeQuantities();

	/**
	 * Compute the divergence of these normalized grad u vectors, at vertex locations
	 * @note can cache edge normal cross product for first loop, and edges and cotangents for second loop
	 */
	void updateDivx();

	/**
	 * Do a pass to compute face gradients of u with the current solution
	 */
	void updateGradu();

	/**
	 * Solves the distances Uses Poisson equation,
	 * Laplacian of distance equal to divergence of normalized heat gradients.
	 * This is step III in Algorithm 1 ofthe Geodesics in Heat paper,
	 * but here is done iteratively with a Gauss-Seidelsolve of some number of steps
	 * to refine the solution whenever this method is called.
	 * @param number of Gauss-Seidel steps to take
	 */
	void solveDistanceStep(int GSSteps);

	/**
	 * Computes the cotangent Laplacian weights at each vertex
	 * Assume no boundaries for now...
	 * Also must assume a triangular mesh!
	 */
	void computeLaplacian();

	/** computes the angle between this he and the next */
	double angleWithNext(HalfEdge *he);
};
#endif