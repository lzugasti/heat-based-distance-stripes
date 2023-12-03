// LEOPOLDO ZUGASTI 260919951

#include "HEDS.h"

HEDS::HEDS(shared_ptr<PolygonSoup> soup)
{
    halfEdges->clear();
    faces->clear();
    faces->reserve(soup->faceList.size());
    vertices = soup->vertexList;
    for (auto &face : soup->faceList) {
		/**
		 * 2 Build the half edge data structure from the polygon soup, triangulating non-triangular faces.
		 */
        for (int i = 0; i < face.size() - 2; i++)
        {
            HalfEdge* he1 = createHalfEdge(soup, face.at(0), face.at(i + 1));
            HalfEdge* he2 = createHalfEdge(soup, face.at(i + 1), face.at(i + 2));
            HalfEdge* he3 = createHalfEdge(soup, face.at(i + 2), face.at(0));
            he1->next = he2;
            he2->next = he3;
            he3->next = he1;
            faces->push_back(new Face(he1));
        }
    }
    // set vertex normals
    for (auto &v : *vertices) {
        v->n = glm::vec3(0, 0, 0);
    }
    for (auto &f : *faces) {
        f->computeNormal();
        HalfEdge *he = f->he;
        do
        {
            he->e = he->twin->head->p - he->head->p;
            he->twin->e = -he->e;

            he->ecn = glm::cross(f->n, he->next->next->e);
            he->twin->ecn = glm::cross(he->twin->leftFace->n, he->twin->next->next->e);
            he = he->next;
        } while (he != f->he);
        
    }
    /**
     * 3 Compute vertex normals.
     */
    for (auto& v : *vertices)
    {
        float count = 0;
        HalfEdge* he = v->he;
        do {
            count++;
            v->n += he->leftFace->n;
            he = he->next->twin;
        } while (he != v->he);
        v->n /= glm::length(v->n);
    }
}

HalfEdge *HEDS::createHalfEdge(shared_ptr<PolygonSoup> soup, int i, int j)
{
    std::string p = to_string(i) + "," + to_string(j);
    if (halfEdges->count(p) > 0) {
        throw runtime_error("non orientable manifold");
    }
    std::string twin = to_string(j) + "," + to_string(i);
    HalfEdge *he = new HalfEdge();
    he->head = soup->vertexList->at(j);
    he->head->he = he; // make sure the vertex has at least one half edge that points to it.
    int twinCount = halfEdges->count(twin);
    if (twinCount > 0) {
        he->twin = halfEdges->at(twin);
        he->twin->twin = he;
    }
    halfEdges->emplace(p, he);
    return he;
}

void HEDS::initHeatFlow()
{
    for (auto &v : *vertices) {
        if (v->constrained) {
            v->u0 = 1;
            v->ut = 1;
        } else {
            v->u0 = 0;
            v->ut = 0;
        }
        v->phi = 0;
    }
}

void HEDS::solveHeatFlowStep(int GSSteps, double t)
{
    // we'll naively choose some random vertex as a source, and
    // then do lots of GS iterations on that for a backward Euler solve of
    // (A-tL) u_t = u_0
    // what is a good value of t?  t small for accuracy, but t large for possibly floating point error
    for (int i = 0; i < GSSteps; i++) {
        for (auto &v : *vertices) {
            if (v->constrained) {
                continue; // do nothing for the constrained vertex!
            }
            /**
             * 7 write inner loop code for the PGS heat solve.
             */
            double sum = 0;
            int j = 0;
            HalfEdge* h = v->he;
            do {
                sum += (h->twin->head->area - t*v->Lij[j]) * h->twin->head->ut;
                h = h->next->twin;
                j++;
            } while (h != v->he);
            double val = (1.0 / (v->area - t*v->Lii)) * (v->u0 - sum);
            v->ut = val;
            
        }
    }
}

void HEDS::precomputeQuantities()
{
    /**
     * TODO: you can do some pre-computation here to make things faster!
     */

    for (auto &f : *faces) {
        HalfEdge *he = f->he;
    }

    for (auto &v : *vertices) {
        v->divX = 0;
    }
}

void HEDS::updateDivx()
{
    // Compute the divergence of these normalized grad u vectors, at vertex locations
    for (auto &v : *vertices) {
        v->divX = 0;
        /**
         * TODO: 9 Update the divergence of the normalized gradients, ie., v.divX for each Vertex v.
         */
        HalfEdge* h = v->he;
        do
        {
            float cotan1 = 1/glm::tan(angleWithNext(h->next->next));
            float cotan2 = 1/glm::tan(angleWithNext(h->next));
            glm::vec3 e1 = h->twin->e;
            glm::vec3 e2 = h->next->e;
            glm::vec3 X = h->leftFace->gradu;
            if (isnan(X[0]) || isinf(X[0])) {
                h = h->next->twin;
                continue;
            }
            v->divX += cotan2 *(glm::dot(e1, X)) + cotan1 *(glm::dot(e2, X));
            h = h->next->twin;
        } while (h != v->he);
        v->divX = v->divX / 2.0;
    }
}

void HEDS::updateGradu() {
    // do a pass to compute face gradients of u with the current solution
    for (auto& f : *faces) {
        f->gradu[0] = 0;
        f->gradu[1] = 0;
        f->gradu[2] = 0;
        /**
         * 8 update the gradient of u from the heat values, i.e., f.gradu for each Face f.
         */
        HalfEdge* h = f->he;
        do
        {
            f->gradu += h->ecn * (float)h->head->ut;
            h = h->next;
        } while (h != f->he);

        f->gradu = glm::normalize(f->gradu);
    }
}

void HEDS::solveDistanceStep(int GSSteps)
{
    // Finally step the solution to the distance problem
    for (int i = 0; i < GSSteps; i++) {
        for (auto &v : *vertices) {
            // LHS matrix is L, so to take all the LHS to the RSH for one variable we get
            // Lii phi_i = div X + sum_{j!=i} tLij phi_j
            /**
             * TODO: 10 Implement the inner loop of the Gauss-Seidel solve to compute the distances to each vertex, phi.
             */
            double sum = 0;
            int j = 0;
            HalfEdge* h = v->he;
            do {
                sum += v->Lij[j] * h->twin->head->phi;
                h = h->next->twin;
                j++;
            } while (h != v->he);
            double val = (1 / v->Lii) * (v->divX - sum);
            v->phi = val; 
        }
    }

    // Note that the solution to step III is unique only up to an additive constant,
    // final values simply need to be shifted such that the smallest distance is zero. 
    // We also identify the max phi value here to identify the maximum geodesic and to 
    // use adjusting the colour map for rendering
    minphi = DBL_MAX;
    maxphi = DBL_MIN;
    for (auto &v : *vertices) {
        if (v->phi < minphi)
            minphi = v->phi;
        if (v->phi > maxphi)
            maxphi = v->phi;
    }
    maxphi -= minphi;
    for (auto& v : *vertices) {
        v->phi -= minphi;
    }
}

void HEDS::computeLaplacian()
{
    for (auto &v : *vertices) {
        v->area = 0;
        v->Lii = 0;
        // get degree of v
        int degree = v->valence();
        v->Lij = new double[degree];

        /**
         * 6 Compute the Laplacian and store as vertex weights, and cotan operator diagonal Lii and off diagonal Lij terms. 
         */


        double sum = 0;
        HalfEdge* h = v->he;
        int j = 0;
        do {
            double cotBeta = 1 / glm::tan(angleWithNext(h->next));
            double cotAlpha = 1 / glm::tan(angleWithNext(h->twin->next));

            v->Lij[j] = (cotBeta + cotAlpha) * (- 1 / 2.0f);
            sum += cotBeta + cotAlpha;

            v->area += h->leftFace->area;

            j++;
            h = h->next->twin;
        } while (h != v->he);
        v->Lii = sum * (1 / 2.0f);
        v->area = v->area / 3.0;
    }
}

double HEDS::angleWithNext(HalfEdge *he)
{
    /**
     * 6 Implement this function to compute the angle with next edge... you'll want to use this in a few places.
     */
    return glm::acos(glm::dot(he->twin->e, he->next->twin->e) / (glm::length(he->twin->e) * glm::length(he->next->twin->e)));
}
