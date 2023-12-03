// LEOPOLDO ZUGASTI 260919951

#pragma once
#ifndef Face_H
#define Face_H
#include <memory>

#include "GLHeaders.h"

#include "HalfEdge.h"

using namespace std;

/**
 * Simple face class
 */
class Face
{
public:
    /** sure, why not keep a normal for flat shading? */
    glm::vec3 n;

    /** Some half edge on the face */
    HalfEdge *he;

    /** area of the face */
    double area;

    /** gradient of the quantities u at the vertices */
    glm::vec3 gradu;

    /** center of the face */
    glm::vec3 c;

    /**
     * Constructs a face from a half edge, and computes the flat normal
     * @param he
     */
    Face(HalfEdge *he);

    /**
     * Computes the center of the face
     */
    void computeCenter();

    /**
     * Computes the non-normalized normal vector
     */
    void computeNormal();
};
#endif