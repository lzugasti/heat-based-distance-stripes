// LEOPOLDO ZUGASTI 260919951

#include "MeshDrawHeatGeo.h"
void MeshDrawHeatGeo::init(shared_ptr<Program> pProgram)
{
    this->aProgram = pProgram;

    camSpacePositionID = aProgram->getAttribute("vertex");
    camSpaceNormalID = aProgram->getAttribute("normal");
    utvID = aProgram->getAttribute("ut");
    phivID = aProgram->getAttribute("phi");
}

void MeshDrawHeatGeo::setMesh(shared_ptr<HEDS> heds)
{
    this->heds = heds;

    glGenVertexArrays(1, &vao);

    pvdCount = (int)heds->vertices->size();
    int pvdSize = 3 + 3; // vertex, normal
    int dvdSize = 1 + 1; // u and phi
    pvdStrideBytes = pvdSize * sizeof(GLfloat);
    dvdStrideBytes = dvdSize * sizeof(GLfloat);

    if (pvdBuffer != NULL)
    {
        free(pvdBuffer);
    }
    if (dvdBuffer != NULL)
    {
        free(dvdBuffer);
    }
    if (indexBuffer != NULL)
    {
        free(indexBuffer);
    }

    pvdBuffer = new GLfloat[pvdSize * pvdCount];
    dvdBuffer = new GLfloat[dvdSize * pvdCount];
    indexBuffer = new GLuint[heds->faces->size() * 3];

    fillStaticVBOs();
    fillDynamicVBO();

    glGenBuffers(1, &pvdVBOID);
    glBindBuffer(GL_ARRAY_BUFFER, pvdVBOID);
    glBufferData(GL_ARRAY_BUFFER, pvdStrideBytes * pvdCount, pvdBuffer, GL_STATIC_DRAW);

    glGenBuffers(1, &dvdVBOID);
    glBindBuffer(GL_ARRAY_BUFFER, dvdVBOID);
    glBufferData(GL_ARRAY_BUFFER, dvdStrideBytes * pvdCount, dvdBuffer, GL_DYNAMIC_DRAW);

    glGenBuffers(1, &indexVBOID);
    glBindBuffer(GL_ARRAY_BUFFER, indexVBOID);
    glBufferData(GL_ARRAY_BUFFER, heds->faces->size() * 3 * sizeof(GLuint), indexBuffer, GL_STATIC_DRAW);

    // Unbind the arrays
    glBindBuffer(GL_ARRAY_BUFFER, 0);

    if (gradientBuffer != NULL)
    {
        free(gradientBuffer);
    }

    gradientBufferSize = heds->faces->size() * 3 * 2;
    gradientBuffer = new GLfloat[gradientBufferSize];

    glGenVertexArrays(1, &gradientVao);

    glGenBuffers(1, &gradientBufferID);
    glBindBuffer(GL_ARRAY_BUFFER, gradientBufferID);
    glBufferData(GL_ARRAY_BUFFER, gradientBufferSize * sizeof(GLfloat), gradientBuffer, GL_DYNAMIC_DRAW);

    // Unbind the arrays
    glBindBuffer(GL_ARRAY_BUFFER, 0);
}

void MeshDrawHeatGeo::fillDynamicVBO()
{
    double minu = DBL_MAX;
    double maxu = DBL_MIN; // can just be 1
    for (auto &v : *(heds->vertices)) {
        if (v->ut < minu)
            minu = v->ut;
        if (v->ut > maxu)
            maxu = v->ut;
    }

    double s = scale;
    int count = 0;
    for (auto &v : *(heds->vertices)) {
        double val = sqrt(-log(max(v->ut, 1e-13)));
        double valmax = sqrt(-log(max(minu, 1e-13)));
        val = 1 - 1.5 * val / valmax;

        dvdBuffer[count] = (GLfloat)val; // Varadahn's formula
        count++;
        dvdBuffer[count] = (GLfloat)(s * (v->phi - heds->minphi) / (heds->maxphi - heds->minphi));
        count++;
    }
}

void MeshDrawHeatGeo::fillStaticVBOs()
{
    int count = 0;
    for (auto &v : *(heds->vertices))
    {
        pvdBuffer[count] = (GLfloat)v->p.x;
        count++;
        pvdBuffer[count] = (GLfloat)v->p.y;
        count++;
        pvdBuffer[count] = (GLfloat)v->p.z;
        count++;
        pvdBuffer[count] = (GLfloat)v->n.x;
        count++;
        pvdBuffer[count] = (GLfloat)v->n.y;
        count++;
        pvdBuffer[count] = (GLfloat)v->n.z;
        count++;
    }

    count = 0;
    for (auto &f : *(heds->faces))
    {
        indexBuffer[count] = f->he->head->index;
        count++;
        indexBuffer[count] = f->he->next->head->index;
        count++;
        indexBuffer[count] = f->he->next->next->head->index;
        count++;
    }
}

void MeshDrawHeatGeo::drawVBOs(shared_ptr<MatrixStack> P, shared_ptr<MatrixStack> MV)
{
    aProgram->bind();

    glBindVertexArray(vao);

    glUniformMatrix4fv(aProgram->getUniform("P"), 1, GL_FALSE, &P->topMatrix()[0][0]);
    glUniformMatrix4fv(aProgram->getUniform("MV"), 1, GL_FALSE, &MV->topMatrix()[0][0]);

    const GLfloat lightCamSpacePosition[3] = {0, 20, 20};
    const GLfloat lightColor[3] = {1, 1, 1};
    const GLfloat materialDiffuse[3] = {0.2f, 0.2f, 1.0f};

    glUniform3fv(aProgram->getUniform("lightCamSpacePosition"), 1, lightCamSpacePosition);
    glUniform3fv(aProgram->getUniform("lightColor"), 1, lightColor);
    glUniform3fv(aProgram->getUniform("materialDiffuse"), 1, materialDiffuse);
    glUniform1f(aProgram->getUniform("materialShininess"), 255);

    glEnableVertexAttribArray(camSpacePositionID);
    glEnableVertexAttribArray(camSpaceNormalID);
    glEnableVertexAttribArray(utvID);
    glEnableVertexAttribArray(phivID);

    glBindBuffer(GL_ARRAY_BUFFER, pvdVBOID);
    glVertexAttribPointer(camSpacePositionID, 3, GL_FLOAT, false, pvdStrideBytes, (void *)(0 * 4)); // last parameter is byte offset
    glVertexAttribPointer(camSpaceNormalID, 3, GL_FLOAT, false, pvdStrideBytes, (void *)(3 * 4));

    fillDynamicVBO(); // refill the vertex data for heat diffusion info
    glBindBuffer(GL_ARRAY_BUFFER, dvdVBOID);
    glBufferSubData(GL_ARRAY_BUFFER, 0, dvdStrideBytes * pvdCount, dvdBuffer); // use substitution
    glVertexAttribPointer(utvID, 1, GL_FLOAT, false, dvdStrideBytes, (void *)(0 * 4)); // last parameter is byte offset
    glVertexAttribPointer(phivID, 1, GL_FLOAT, false, dvdStrideBytes, (void *)(1 * 4));

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, indexVBOID);
    glDrawElements(GL_TRIANGLES, heds->faces->size() * 3, GL_UNSIGNED_INT, 0);

    glDisableVertexAttribArray(camSpacePositionID);
    glDisableVertexAttribArray(camSpaceNormalID);
    glDisableVertexAttribArray(utvID);
    glDisableVertexAttribArray(phivID);

    glBindBuffer(GL_ARRAY_BUFFER, 0);

    aProgram->unbind();
}

void MeshDrawHeatGeo::drawGrad(shared_ptr<MatrixStack> P, shared_ptr<MatrixStack> MV, shared_ptr<Program> lineProg)
{
    lineProg->bind();

    glBindVertexArray(gradientVao);

    glUniformMatrix4fv(lineProg->getUniform("P"), 1, GL_FALSE, &P->topMatrix()[0][0]);
    glUniformMatrix4fv(lineProg->getUniform("MV"), 1, GL_FALSE, &MV->topMatrix()[0][0]);

    glUniform3f(lineProg->getUniform("col"), 1, 0, 0);

    // Bind position buffer
    int h_pos = lineProg->getAttribute("aPos");
    glEnableVertexAttribArray(h_pos);
    glBindBuffer(GL_ARRAY_BUFFER, gradientBufferID);
    glVertexAttribPointer(h_pos, 3, GL_FLOAT, GL_FALSE, 0, (const void *)0);

    float h = gradHeight;
    float s = gradScale;

    int bufferIndexCount = 0;
    for (auto &f : *(heds->faces))
    {
        gradientBuffer[bufferIndexCount] = f->c.x + f->n.x * h;
        gradientBuffer[bufferIndexCount + 1] = f->c.y + f->n.y * h;
        gradientBuffer[bufferIndexCount + 2] = f->c.z + f->n.z * h;
        gradientBuffer[bufferIndexCount + 3] = f->c.x + f->n.x * h + f->gradu.x * s;
        gradientBuffer[bufferIndexCount + 4] = f->c.y + f->n.y * h + f->gradu.y * s;
        gradientBuffer[bufferIndexCount + 5] = f->c.z + f->n.z * h + f->gradu.z * s;
        bufferIndexCount += 6;
    }

    glBufferSubData(GL_ARRAY_BUFFER, 0, gradientBufferSize * sizeof(GLfloat), &gradientBuffer[0]);

    int pointCount = 0;
    for (auto& f : *(heds->faces))
    {
        glDrawArrays(GL_LINES, pointCount, 2);
        pointCount += 2;
    }

    glDisableVertexAttribArray(h_pos);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    lineProg->unbind();
}
