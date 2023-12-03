// LEOPOLDO ZUGASTI 260919951

#include "GLHeaders.h"

#include "HalfEdge.h"

int Vertex::valence() {
    int v = 0;

    /**
     * 5 compute the valence of this vertex
     */

    HalfEdge* h = he;
    do {
        v++;
        h = h->next->twin;
    } while (h != he);

    return v;
}