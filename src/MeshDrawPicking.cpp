// LEOPOLDO ZUGASTI 260919951

#include "MeshDrawPicking.h"

void MeshDrawPicking::init(shared_ptr<Program> pProgram)
{
    this->aProgram = pProgram;
    this->vertexID = pProgram->getAttribute("vertex");
    this->colorID = pProgram->getAttribute("color");
    pickedTriangleVerticesBuffer.clear();
    pickedTriangleVerticesBuffer.resize(9);
}

void MeshDrawPicking::setMesh(shared_ptr<HEDS> heds)
{
    this->heds = heds;
    nFaces = (int)heds->faces->size();

    pickingVertexB.clear();
    pickingVertexB.reserve(3 * nFaces * 3);
    pickingFaceIDBuffer.clear();
    pickingFaceIDBuffer.reserve(3 * nFaces * 3);
    pickingTriangleIndexBuffer.clear();
    pickingTriangleIndexBuffer.reserve(3 * nFaces);
    fillPickingBuffers();

    glGenVertexArrays(1, &pickVao);

    glGenBuffers(1, &pickingVertexBOID);
    glBindBuffer(GL_ARRAY_BUFFER, pickingVertexBOID);
    glBufferData(GL_ARRAY_BUFFER, pickingVertexB.size() * sizeof(GLfloat), pickingVertexB.data(), GL_STATIC_DRAW);

    glGenBuffers(1, &pickingFaceIDBOID);
    glBindBuffer(GL_ARRAY_BUFFER, pickingFaceIDBOID);
    glBufferData(GL_ARRAY_BUFFER, pickingFaceIDBuffer.size() * sizeof(GLubyte), pickingFaceIDBuffer.data(), GL_STATIC_DRAW);

    glGenBuffers(1, &pickingTriangleIndexBOID);
    glBindBuffer(GL_ARRAY_BUFFER, pickingTriangleIndexBOID);
    glBufferData(GL_ARRAY_BUFFER, pickingTriangleIndexBuffer.size() * sizeof(GLuint), pickingTriangleIndexBuffer.data(), GL_STATIC_DRAW);

    glGenBuffers(1, &pickedVBOID);
    glBindBuffer(GL_ARRAY_BUFFER, pickedVBOID);
    glBufferData(GL_ARRAY_BUFFER, pickedTriangleVerticesBuffer.size() * sizeof(GLfloat), pickedTriangleVerticesBuffer.data(), GL_DYNAMIC_DRAW);
}

void MeshDrawPicking::fillPickingBuffers()
{
    int faceIndex = 0;
    for (auto &f : *(heds->faces))
    {
        GLubyte r = (GLubyte)(faceIndex & 0xff);
        GLubyte g = (GLubyte)((faceIndex >> 8) & 0xff);
        GLubyte b = (GLubyte)((faceIndex >> 16) & 0xff);
        HalfEdge *he = f->he;
        int i = 0;
        do
        {
            Vertex *v = he->head;

            pickingVertexB.push_back((GLfloat)v->p.x);
            pickingVertexB.push_back((GLfloat)v->p.y);
            pickingVertexB.push_back((GLfloat)v->p.z);
            pickingFaceIDBuffer.push_back(r);
            pickingFaceIDBuffer.push_back(g);
            pickingFaceIDBuffer.push_back(b);

            pickingTriangleIndexBuffer.push_back((GLuint)(faceIndex * 3 + i));
            i++;
            he = he->next;
        } while (he != f->he);
        faceIndex++;
    }
    faceIndex = faceIndex;
}

int MeshDrawPicking::pickTriangles(glm::vec2 *p, glm::vec3 *BCC, GLFWwindow *window, shared_ptr<MatrixStack> P, shared_ptr<MatrixStack> MV, shared_ptr<Program> triProg)
{
    int ID = -1;

    glClearColor(1, 1, 1, 1);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    aProgram->bind();

    glBindVertexArray(pickVao);

    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

    glUniformMatrix4fv(aProgram->getUniform("P"), 1, GL_FALSE, &P->topMatrix()[0][0]);
    glUniformMatrix4fv(aProgram->getUniform("MV"), 1, GL_FALSE, &MV->topMatrix()[0][0]);

    glEnableVertexAttribArray(vertexID);
    glEnableVertexAttribArray(colorID);

    glBindBuffer(GL_ARRAY_BUFFER, pickingVertexBOID);
    glVertexAttribPointer(vertexID, 3, GL_FLOAT, false, pickingVertexBOStrideBytes, 0);
    glBindBuffer(GL_ARRAY_BUFFER, pickingFaceIDBOID);
    glVertexAttribPointer(colorID, 3, GL_UNSIGNED_BYTE, false, pickingFaceIDBOStrideBytes, 0);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, pickingTriangleIndexBOID);
    glDrawElements(GL_TRIANGLES, 3 * nFaces, GL_UNSIGNED_INT, 0);

    glDisableVertexAttribArray(vertexID);
    glDisableVertexAttribArray(colorID);

    glBindBuffer(GL_ARRAY_BUFFER, 0);

    aProgram->unbind();

    GLubyte colorPixels[4];

    if (p != NULL)
    {
        // Get current frame buffer size.
        int width, height;
        glfwGetFramebufferSize(window, &width, &height);

        int x = (int)p->x;
        int y = height - (int)p->y;
        glReadPixels(x, y, 1, 1, GL_RGB, GL_UNSIGNED_BYTE, colorPixels);
        int r, g, b;
        r = colorPixels[0] & 0xff;
        g = colorPixels[1] & 0xff;
        b = colorPixels[2] & 0xff;
        ID = r + (g << 8) + (b << 16);

        if (ID >= nFaces)
        {
            return -1;
        }

        triProg->bind();

        glDisable(GL_DEPTH_TEST);

        glUniformMatrix4fv(triProg->getUniform("P"), 1, GL_FALSE, &P->topMatrix()[0][0]);
        glUniformMatrix4fv(triProg->getUniform("MV"), 1, GL_FALSE, &MV->topMatrix()[0][0]);

        int h_pos = triProg->getAttribute("aPos");
        glEnableVertexAttribArray(h_pos);
        glBindBuffer(GL_ARRAY_BUFFER, pickedVBOID);
        glVertexAttribPointer(h_pos, 3, GL_FLOAT, GL_FALSE, 0, (const void *)0);

        glUniform3f(triProg->getUniform("col"), 0, 1, 0);

        HalfEdge *he = heds->faces->at(ID)->he;
        glm::vec3 v1 = he->head->p;
        glm::vec3 v2 = he->next->head->p;
        glm::vec3 v3 = he->next->next->head->p;

        pickedTriangleVerticesBuffer[0] = v1.x;
        pickedTriangleVerticesBuffer[1] = v1.y;
        pickedTriangleVerticesBuffer[2] = v1.z;
        pickedTriangleVerticesBuffer[3] = v2.x;
        pickedTriangleVerticesBuffer[4] = v2.y;
        pickedTriangleVerticesBuffer[5] = v2.z;
        pickedTriangleVerticesBuffer[6] = v3.x;
        pickedTriangleVerticesBuffer[7] = v3.y;
        pickedTriangleVerticesBuffer[8] = v3.z;

        glBufferSubData(GL_ARRAY_BUFFER, 0, 9 * sizeof(GLfloat), pickedTriangleVerticesBuffer.data());

        glDrawArrays(GL_TRIANGLES, 0, 3);

        glDisableVertexAttribArray(h_pos);

        glBindBuffer(GL_ARRAY_BUFFER, 0);

        glEnable(GL_DEPTH_TEST);

        triProg->unbind();

        r = colorPixels[0] & 0xff;
        g = colorPixels[1] & 0xff;
        b = colorPixels[2] & 0xff;
        BCC->x = (float)r / 255.0f;
        BCC->y = (float)g / 255.0f;
        BCC->z = (float)b / 255.0f;
    }
    return ID;
}
