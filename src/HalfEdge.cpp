// LEOPOLDO ZUGASTI 260919951

#include "HalfEdge.h"

HalfEdge *HalfEdge::prev()
{
    HalfEdge *prev = this;
    while (prev->next != this)
        prev = prev->next;
    return prev;
}

GLfloat halfEdgeVerticies[3][3];

void HalfEdge::display(shared_ptr<Program> prog, shared_ptr<MatrixStack> P, shared_ptr<MatrixStack> MV, GLuint vao, GLuint vbo)
{

    //// Draw lines.
    prog->bind();

    glUniformMatrix4fv(prog->getUniform("P"), 1, GL_FALSE, &P->topMatrix()[0][0]);
    glUniformMatrix4fv(prog->getUniform("MV"), 1, GL_FALSE, &MV->topMatrix()[0][0]);
    glBindVertexArray(vao);
    glUniform3f(prog->getUniform("col"), 0, 0.8f, 0);

    glm::vec3 p0 = prev()->head->p;
    glm::vec3 p1 = head->p;
    glm::vec3 p2 = next->head->p;
    double x, y, z;

    x = p0.x * 0.8 + (p1.x + p2.x) * 0.1;
    y = p0.y * 0.8 + (p1.y + p2.y) * 0.1;
    z = p0.z * 0.8 + (p1.z + p2.z) * 0.1;

    halfEdgeVerticies[0][0] = x;
    halfEdgeVerticies[0][1] = y;
    halfEdgeVerticies[0][2] = z;

    x = p1.x * 0.8 + (p0.x + p2.x) * 0.1;
    y = p1.y * 0.8 + (p0.y + p2.y) * 0.1;
    z = p1.z * 0.8 + (p0.z + p2.z) * 0.1;

    halfEdgeVerticies[1][0] = x;
    halfEdgeVerticies[1][1] = y;
    halfEdgeVerticies[1][2] = z;

    x = p1.x * 0.7 + p0.x * 0.1 + p2.x * 0.2;
    y = p1.y * 0.7 + p0.y * 0.1 + p2.y * 0.2;
    z = p1.z * 0.7 + p0.z * 0.1 + p2.z * 0.2;

    halfEdgeVerticies[2][0] = x;
    halfEdgeVerticies[2][1] = y;
    halfEdgeVerticies[2][2] = z;

    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    glLineWidth(1);
    glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(halfEdgeVerticies), &halfEdgeVerticies);
    glUniform3f(prog->getUniform("col"), 0, 1, 0);
    glDrawArrays(GL_LINES, 0, 2);
    glUniform3f(prog->getUniform("col"), 0, 1, 0);
    glDrawArrays(GL_LINES, 1, 2);
    prog->unbind();
}