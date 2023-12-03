// LEOPOLDO ZUGASTI 260919951

#include "Face.h"

using namespace std;

Face::Face(HalfEdge *he)
{
    this->he = he;
    HalfEdge *loop = he;
    do
    {
        loop->leftFace = this;
        loop = loop->next;
    } while (loop != he);
    computeNormal();
    computeCenter();
}

void Face::computeCenter()
{
    c = glm::vec3(0, 0, 0);
    HalfEdge *loop = he;
    do
    {
        c += loop->head->p;
        loop = loop->next;
    } while (loop != he);
    c *= (1.0 / 3.0);
}

void Face::computeNormal()
{
    glm::vec3 p0 = he->head->p;
    glm::vec3 p1 = he->next->head->p;
    glm::vec3 p2 = he->next->next->head->p;
    glm::vec3 v1;
    glm::vec3 v2;
    v1 = p1 - p0;
    v2 = p2 - p1;
    n = glm::cross(v1, v2);
    //area = 0.5 * glm::length(n);
    area = 0.5 * n.length();
    n = glm::normalize(n);
}