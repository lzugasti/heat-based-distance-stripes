// LEOPOLDO ZUGASTI 260919951

#pragma once
#ifndef PolygonSoup_H
#define PolygonSoup_H

#include <vector>
#include <map>
#include <string>

#include <iostream>
#include <sstream>
#include <fstream>

#include <float.h>

#include "GLHeaders.h"
#include "GLSL.h"
#include "Program.h"
#include "Vertex.h"
#include "MatrixStack.h"

using namespace std;

/**
 * Simple implementation of a loader for a polygon soup
 */
class PolygonSoup
{
public:
    ~PolygonSoup()
    {
        posBuf.clear();
        norBuf.clear();
        faceList.clear();
        faceSidesHistogram.clear();
        vertexList->clear();
        vertexList = NULL;
        prog = NULL;
    };

    /* position buffer */
    std::vector<float> posBuf;
    /* normal buffer */
    std::vector<float> norBuf;
    /* local vao */
    GLuint vao;
    /* position buffer ID */
    GLuint posBufID;
    /* normal buffer ID */
    GLuint norBufID;

    /**
     * List of vertex objects used in the mes
     */
    shared_ptr<vector<Vertex *>> vertexList = make_shared<vector<Vertex *>>();

    /**
     * List of faces, where each face is a list indices into the vertex list.
     */
    vector<vector<int>> faceList;

    /* default drawing program */
    shared_ptr<Program> prog;

    /* soup face sides histogram */
    map<int, int> faceSidesHistogram;

    /* soup statistics string */
    string soupStatistics;

    /**
     * Creates a polygon soup by loading an OBJ file
     * @param file
     */
    PolygonSoup(string filename);

    /**
     * Parses a vertex definition from a line in an obj file, and
     * directly inserts it into the vertex list.
     * Assumes that there are three components.
     * @param newline
     */
    Vertex *parseVertex(string newline);

    /**
     * Gets the list of indices for a face from a string in an obj file.
     * Simply ignores texture and normal information for simplicity
     * @param newline
     * @return list of indices
     */

    vector<int> parseFace(string newline);

    /**
     * Draw the polygon soup
     */
    void display(shared_ptr<Program> prog, shared_ptr<MatrixStack> P, shared_ptr<MatrixStack> MV);

    /**
     * Get a line safely
     * @see http://stackoverflow.com/questions/6089231/getting-std-ifstream-to-handle-lf-cr-and-crlf
     */
    static istream &safeGetline(istream &is, string &t)
    {
        t.clear();

        // The characters in the stream are read one-by-one using a streambuf.
        // That is faster than reading them one-by-one using the istream.
        // Code that uses streambuf this way must be guarded by a sentry object.
        // The sentry object performs various tasks,
        // such as thread synchronization and updating the stream state.

        istream::sentry se(is, true);
        streambuf *sb = is.rdbuf();

        if (se)
        {
            for (;;)
            {
                int c = sb->sbumpc();
                switch (c)
                {
                case '\n':
                    return is;
                case '\r':
                    if (sb->sgetc() == '\n')
                        sb->sbumpc();
                    return is;
                case EOF:
                    // Also handle the case when the last line has no line ending
                    if (t.empty())
                        is.setstate(ios::eofbit);
                    return is;
                default:
                    t += static_cast<char>(c);
                }
            }
        }

        return is;
    }
};
#endif