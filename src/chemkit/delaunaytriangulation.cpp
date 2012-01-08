/******************************************************************************
**
** Copyright (C) 2009-2011 Kyle Lutz <kyle.r.lutz@gmail.com>
** All rights reserved.
**
** This file is a part of the chemkit project. For more information
** see <http://www.chemkit.org>.
**
** Redistribution and use in source and binary forms, with or without
** modification, are permitted provided that the following conditions
** are met:
**
**   * Redistributions of source code must retain the above copyright
**     notice, this list of conditions and the following disclaimer.
**   * Redistributions in binary form must reproduce the above copyright
**     notice, this list of conditions and the following disclaimer in the
**     documentation and/or other materials provided with the distribution.
**   * Neither the name of the chemkit project nor the names of its
**     contributors may be used to endorse or promote products derived
**     from this software without specific prior written permission.
**
** THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
** "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
** LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
** A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
** OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
** SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
** LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
** DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
** THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
** (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
** OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
**
******************************************************************************/

#include "delaunaytriangulation.h"

#include <set>
#include <deque>
#include <vector>
#include <algorithm>

#include "point3.h"
#include "foreach.h"
#include "vector3.h"
#include "geometry.h"
#include "alphashape.h"

namespace chemkit {

namespace {

// === EdgeSet ============================================================= //
class EdgeSet
{
public:
    EdgeSet(int vertexCount);

    void insert(int a, int b);
    bool contains(int a, int b);

private:
    std::vector<std::set<int> > m_edges;
};

EdgeSet::EdgeSet(int vertexCount)
    : m_edges(vertexCount)
{
}

void EdgeSet::insert(int a, int b)
{
    if(a > b)
        std::swap(a, b);

    m_edges[a].insert(b);
}

bool EdgeSet::contains(int a, int b)
{
    if(a > b)
        std::swap(a, b);

    return m_edges[a].count(b) != 0;
}

// === Tetrahedron ========================================================= //
class Tetrahedron
{
public:
    int verticies[4];
    int neighbors[4];
    bool valid;
    bool inAlphaShape;

    bool contains(int vertex) const;
    DelaunayTriangulation::Triangle triangle(int index) const;
};

bool Tetrahedron::contains(int vertex) const
{
    for(int i = 0; i < 4; i++){
        if(verticies[i] == vertex){
            return true;
        }
    }

    return false;
}

DelaunayTriangulation::Triangle Tetrahedron::triangle(int index) const
{
    DelaunayTriangulation::Triangle triangle;

    // abc
    if(index == 0){
        triangle[0] = verticies[0];
        triangle[1] = verticies[1];
        triangle[2] = verticies[2];
    }
    // adb
    else if(index == 1){
        triangle[0] = verticies[0];
        triangle[1] = verticies[3];
        triangle[2] = verticies[1];
    }
    // acd
    else if(index == 2){
        triangle[0] = verticies[0];
        triangle[1] = verticies[2];
        triangle[2] = verticies[3];
    }
    // bdc
    else if(index == 3){
        triangle[0] = verticies[1];
        triangle[1] = verticies[3];
        triangle[2] = verticies[2];
    }

    return triangle;
}

} // end anonymous namespace

// === DelaunayTriangulationPrivate ======================================== //
class DelaunayTriangulationPrivate
{
public:
    std::vector<Point3> verticies;
    std::vector<Real> weights;
    std::vector<Tetrahedron> tetrahedra;

    bool alphaShapeCalculated;

    std::vector<DelaunayTriangulation::Edge> delaunayEdges;
    std::vector<DelaunayTriangulation::Triangle> delaunayTriangles;
    std::vector<std::vector<int> > delaunayTetrahedra;

    std::vector<DelaunayTriangulation::Edge> alphaShapeEdges;
    std::vector<DelaunayTriangulation::Triangle> alphaShapeTriangles;
    std::vector<std::vector<int> > alphaShapeTetrahedra;
};

// === DelaunayTriangulation =============================================== //
/// \class DelaunayTriangulation delaunaytriangulation.h chemkit/delaunaytriangulation.h
/// \ingroup chemkit
/// \internal
/// \brief The DelaunayTriangulation class represents a
///        three-dimensional delaunay triangulation.
///
/// The DelaunayTriangulation class computes and stores the
/// three-dimensional delaunay triangulation.
///
/// The delaunay triangulation is the geometric dual of the
/// voronoi diagram.

// --- Construction and Destruction ---------------------------------------- //
/// Creates a new delaunay triangulation for \p points.
DelaunayTriangulation::DelaunayTriangulation(const std::vector<Point3> &points)
    : d(new DelaunayTriangulationPrivate)
{
    d->verticies = points;

    d->alphaShapeCalculated = false;

    triangulate(false);
}

/// Creates a new weighted delaunay triangulation for \p points with
/// \p weights.
DelaunayTriangulation::DelaunayTriangulation(const std::vector<Point3> &points, const std::vector<Real> &weights)
    : d(new DelaunayTriangulationPrivate)
{
    d->verticies = points;
    d->weights = weights;

    d->alphaShapeCalculated = false;

    triangulate(true);
}

/// Destroys the delaunay triangulation object.
DelaunayTriangulation::~DelaunayTriangulation()
{
    delete d;
}

// --- Properties ---------------------------------------------------------- //
/// Returns the number of points in the delaunay triangulation.
int DelaunayTriangulation::size() const
{
    return vertexCount();
}

/// Returns the position of \p vertex.
Point3 DelaunayTriangulation::position(int vertex) const
{
    return d->verticies[vertex];
}

/// Returns the weight of \p vertex.
Real DelaunayTriangulation::weight(int vertex) const
{
    return d->weights[vertex];
}

/// Returns \c true if the delaunay triangulation has weighted
/// verticies.
bool DelaunayTriangulation::isWeighted() const
{
    return !d->weights.empty();
}

// --- Simplicies ---------------------------------------------------------- //
/// Returns a list of verticies in the delaunay triangulation.
std::vector<int> DelaunayTriangulation::verticies() const
{
    std::vector<int> verticies;

    for(unsigned int i = 0; i < d->verticies.size() - 4; i++){
        verticies.push_back(i);
    }

    return verticies;
}

/// Returns the number of verticies in the delaunay triangulation.
int DelaunayTriangulation::vertexCount() const
{
    return verticies().size();
}

/// Returns a list of edges in the delaunay triangulation.
const std::vector<DelaunayTriangulation::Edge>& DelaunayTriangulation::edges() const
{
    if(d->delaunayEdges.empty()){
        std::vector<Edge> edges;
        EdgeSet edgeSet(d->verticies.size());

        foreach(const std::vector<int> &tetrahedron, tetrahedra()){
            for(int i = 0; i < 4; i++){
                for(int j = i + 1; j < 4; j++){
                    Edge edge;
                    edge[0] = tetrahedron[i];
                    edge[1] = tetrahedron[j];

                    if(!edgeSet.contains(edge[0], edge[1])){
                        edges.push_back(edge);
                        edgeSet.insert(edge[0], edge[1]);
                    }
                }
            }
        }

        d->delaunayEdges = edges;
    }

    return d->delaunayEdges;
}

/// Returns the number of edges in the delaunay triangulation.
int DelaunayTriangulation::edgeCount() const
{
    return edges().size();
}

/// Returns a list of faces in the delaunay triangulation.
const std::vector<DelaunayTriangulation::Triangle>& DelaunayTriangulation::triangles() const
{
    if(d->delaunayTriangles.empty()){
        int initialTetrahedron = 0;
        for(unsigned int i = 0; i < d->tetrahedra.size(); i++){
            const Tetrahedron &tetrahedron = d->tetrahedra[i];
            if(!tetrahedron.valid || isExternal(i)){
                continue;
            }

            initialTetrahedron = i;
            break;
        }

        std::set<int> visited;
        std::deque<int> stack;

        stack.push_front(initialTetrahedron);

        while(!stack.empty()){
            int index = stack.front();
            stack.pop_front();
            visited.insert(index);
            const Tetrahedron &tetrahedron = d->tetrahedra[index];

            for(int i = 0; i < 4; i++){
                int neighborIndex = tetrahedron.neighbors[i];
                if(neighborIndex == -1 || visited.count(neighborIndex)){
                    continue;
                }

                d->delaunayTriangles.push_back(tetrahedron.triangle(i));

                stack.push_front(neighborIndex);
            }
        }
    }

    return d->delaunayTriangles;
}

/// Returns the number of faces in the delaunay triangulation.
int DelaunayTriangulation::triangleCount() const
{
    return triangles().size();
}

/// Returns a list of the tetrahedra in the delaunay triangulation.
const std::vector<std::vector<int> >& DelaunayTriangulation::tetrahedra() const
{
    if(d->delaunayTetrahedra.empty()){
        std::vector<std::vector<int> > tetrahedra;

        foreach(const Tetrahedron &tetrahedron, d->tetrahedra){
            if(!tetrahedron.valid){
                continue;
            }

            std::vector<int> verticies(4);
            bool external = false;
            for(int i = 0; i < 4; i++){
                unsigned int vertex = tetrahedron.verticies[i];

                if(vertex >= (d->verticies.size() - 4)){
                    external = true;
                    break;
                }

                verticies[i] = vertex;
            }

            if(external){
                continue;
            }

            tetrahedra.push_back(verticies);
        }

        d->delaunayTetrahedra = tetrahedra;
    }

    return d->delaunayTetrahedra;
}

/// Returns the number of tetrahedra in the delaunay triangulation.
int DelaunayTriangulation::tetrahedronCount() const
{
    return tetrahedra().size();
}

// --- Geometry ------------------------------------------------------------ //
/// Returns the total volume of the triangulation.
Real DelaunayTriangulation::volume() const
{
    Real volume = 0;

    foreach(const std::vector<int> &tetrahedron, tetrahedra()){
        const Point3 &a = position(tetrahedron[0]);
        const Point3 &b = position(tetrahedron[1]);
        const Point3 &c = position(tetrahedron[2]);
        const Point3 &d = position(tetrahedron[3]);

        volume += chemkit::geometry::tetrahedronVolume(a, b, c, d);
    }

    return volume;
}

/// Returns the total surface area of the triangulation.
Real DelaunayTriangulation::surfaceArea() const
{
    return 0;
}

// --- Alpha Shape --------------------------------------------------------- //
const std::vector<DelaunayTriangulation::Edge>& DelaunayTriangulation::alphaShapeEdges(const AlphaShape *alphaShape) const
{
    if(d->alphaShapeEdges.empty()){
        std::vector<Edge> alphaShapeEdges;
        EdgeSet alphaEdgeSet(d->verticies.size());

        foreach(const Triangle &triangle, alphaShapeTriangles(alphaShape)){
            for(int i = 0; i < 3; i++){
                for(int j = i + 1; j < 3; j++){
                    Edge edge;
                    edge[0] = triangle[i];
                    edge[1] = triangle[j];

                    if(!alphaEdgeSet.contains(edge[0], edge[1])){
                        alphaShapeEdges.push_back(edge);
                        alphaEdgeSet.insert(edge[0], edge[1]);
                    }

                }
            }
        }

        EdgeSet attachedEdgeSet(d->verticies.size());

        foreach(const Triangle &triangle, triangles()){
            int a = triangle[0];
            int b = triangle[1];
            int c = triangle[2];

            if(alphaShape->edgeAttached(a, b, c))
                attachedEdgeSet.insert(a, b);

            if(alphaShape->edgeAttached(a, c, b))
                attachedEdgeSet.insert(a, c);

            if(alphaShape->edgeAttached(b, c, a))
                attachedEdgeSet.insert(b, c);
        }

        foreach(const Edge &edge, edges()){
            if(alphaEdgeSet.contains(edge[0], edge[1])){
                continue;
            }
            else if(attachedEdgeSet.contains(edge[0], edge[1])){
                continue;
            }

            if(alphaShape->orthoradius(edge[0], edge[1]) < alphaShape->alphaValue()){
                alphaShapeEdges.push_back(edge);
            }
        }

        d->alphaShapeEdges = alphaShapeEdges;
    }

    return d->alphaShapeEdges;
}

const std::vector<DelaunayTriangulation::Triangle>& DelaunayTriangulation::alphaShapeTriangles(const AlphaShape *alphaShape) const
{
    if(d->alphaShapeTriangles.empty()){
        calculateAlphaShape(alphaShape);

        std::vector<Triangle> triangles;

        int initialTetrahedron = 0;
        for(unsigned int i = 0; i < d->tetrahedra.size(); i++){
            const Tetrahedron &tetrahedron = d->tetrahedra[i];
            if(!tetrahedron.valid || isExternal(i)){
                continue;
            }

            initialTetrahedron = i;
            break;
        }

        std::set<int> visited;
        std::deque<int> stack;

        stack.push_front(initialTetrahedron);

        while(!stack.empty()){
            int index = stack.front();
            stack.pop_front();
            visited.insert(index);
            const Tetrahedron &tetrahedron = d->tetrahedra[index];

            for(int triangleIndex = 0; triangleIndex < 4; triangleIndex++){
                int neighborIndex = tetrahedron.neighbors[triangleIndex];
                if(neighborIndex == -1 || visited.count(neighborIndex)){
                    continue;
                }

                const Tetrahedron &neighbor = d->tetrahedra[neighborIndex];
                bool neighborExternal = isExternal(neighborIndex);

                if(!neighborExternal){
                    stack.push_front(neighborIndex);
                }
                else{
                    visited.insert(neighborIndex);
                }

                Triangle triangle = tetrahedron.triangle(triangleIndex);

                if(tetrahedron.inAlphaShape && (!neighborExternal && neighbor.inAlphaShape)){
                    triangles.push_back(triangle);
                }
                else if(tetrahedron.inAlphaShape || (!neighborExternal && neighbor.inAlphaShape)){
                    triangles.push_back(triangle);
                }
                else{
                    std::vector<int> tetrahedronVerticies;
                    std::vector<int> neighborVerticies;
                    for(int i = 0; i < 4; i++){
                        tetrahedronVerticies.push_back(tetrahedron.verticies[i]);
                        neighborVerticies.push_back(neighbor.verticies[i]);
                    }

                    for(int i = 0; i < 3; i++){
                        tetrahedronVerticies.erase(std::remove(tetrahedronVerticies.begin(),
                                                               tetrahedronVerticies.end(),
                                                               triangle[i]));

                        neighborVerticies.erase(std::remove(neighborVerticies.begin(),
                                                            neighborVerticies.end(),
                                                            triangle[i]));
                    }

                    int va = triangle[0];
                    int vb = triangle[1];
                    int vc = triangle[2];
                    int vd = tetrahedronVerticies[0];
                    int ve = neighborVerticies[0];

                    if(alphaShape->triangleAttached(va, vb, vc, vd)){
                        continue;
                    }
                    else if(!neighborExternal && alphaShape->triangleAttached(va, vb, vc, ve)){
                        continue;
                    }

                    if(alphaShape->orthoradius(va, vb, vc) < alphaShape->alphaValue()){
                        triangles.push_back(triangle);
                    }
                }
            }
        }

        d->alphaShapeTriangles = triangles;
    }

    return d->alphaShapeTriangles;
}

const std::vector<std::vector<int> >& DelaunayTriangulation::alphaShapeTetrahedra(const AlphaShape *alphaShape) const
{
    if(d->alphaShapeTetrahedra.empty()){
        calculateAlphaShape(alphaShape);

        for(unsigned int i = 0; i < d->tetrahedra.size(); i++){
            const Tetrahedron &tetrahedron = d->tetrahedra[i];
            if(!tetrahedron.valid || isExternal(i))
                continue;

            if(tetrahedron.inAlphaShape){
                std::vector<int> verticies(4);
                for(int j = 0; j < 4; j++){
                    verticies[j] = tetrahedron.verticies[j];
                }

                d->alphaShapeTetrahedra.push_back(verticies);
            }
        }
    }

    return d->alphaShapeTetrahedra;
}

// Marks the tetrahedra that are in the alpha shape.
void DelaunayTriangulation::calculateAlphaShape(const AlphaShape *alphaShape) const
{
    if(d->alphaShapeCalculated){
        return;
    }

    for(unsigned int i = 0; i < d->tetrahedra.size(); i++){
        Tetrahedron &tetrahedron = d->tetrahedra[i];
        if(!tetrahedron.valid || isExternal(i)){
            continue;
        }

        int a = tetrahedron.verticies[0];
        int b = tetrahedron.verticies[1];
        int c = tetrahedron.verticies[2];
        int d = tetrahedron.verticies[3];

        if(alphaShape->orthoradius(a, b, c, d) < alphaShape->alphaValue()){
            tetrahedron.inAlphaShape = true;
        }
        else{
            tetrahedron.inAlphaShape = false;
        }
    }

    d->alphaShapeCalculated = true;
}

// --- Internal Methods ---------------------------------------------------- //
void DelaunayTriangulation::triangulate(bool weighted)
{
    // size of vertex list
    int size = d->verticies.size();

    // build big tetrahedron which will contain all other points. its
    // vertices will be the last four positions in the vertex vector
    d->verticies.push_back(Point3(0, 1e10, 0));
    d->verticies.push_back(Point3(1e10, -1e10, 1e10));
    d->verticies.push_back(Point3(-1e10, -1e10, 1e10));
    d->verticies.push_back(Point3(0, -1e10, -1e10));

    if(weighted){
        d->weights.push_back(0);
        d->weights.push_back(0);
        d->weights.push_back(0);
        d->weights.push_back(0);
    }

    Tetrahedron big;
    big.verticies[0] = size;
    big.verticies[1] = size + 1;
    big.verticies[2] = size + 2;
    big.verticies[3] = size + 3;
    big.neighbors[0] = -1;
    big.neighbors[1] = -1;
    big.neighbors[2] = -1;
    big.neighbors[3] = -1;
    big.valid = true;
    d->tetrahedra.push_back(big);

    // insert verticies
    for(int i = 0; i < size; i++){
        insertPoint(i);
    }
}

/// Returns the index of the tetrahedron that contains the point.
int DelaunayTriangulation::location(const Point3 &point) const
{
    int tetrahedronIndex = 0;

    // find last valid tetrahedron to start at
    for(int i = d->tetrahedra.size() - 1; i != 0; i++){
        const Tetrahedron &tetrahedron = d->tetrahedra[i];
        if(!tetrahedron.valid)
            continue;

        tetrahedronIndex = i;
        break;
    }

    // walk through the delaunay structure and try to find a
    // tetrahedron that contains the point.
    for(size_t iteration = 0; iteration < d->tetrahedra.size(); iteration++){
        const Tetrahedron &tetrahedron = d->tetrahedra[tetrahedronIndex];
        const Point3 &a = position(tetrahedron.verticies[0]);
        const Point3 &b = position(tetrahedron.verticies[1]);
        const Point3 &c = position(tetrahedron.verticies[2]);
        const Point3 &d = position(tetrahedron.verticies[3]);

        if(chemkit::geometry::planeOrientation(a, b, c, point) > 0){
            tetrahedronIndex = tetrahedron.neighbors[0];
        }
        else if(chemkit::geometry::planeOrientation(a, d, b, point) > 0){
            tetrahedronIndex = tetrahedron.neighbors[1];
        }
        else if(chemkit::geometry::planeOrientation(a, c, d, point) > 0){
            tetrahedronIndex = tetrahedron.neighbors[2];
        }
        else if(chemkit::geometry::planeOrientation(b, d, c, point) > 0){
            tetrahedronIndex = tetrahedron.neighbors[3];
        }
        else{
            // we found the tetrahedron that contains the point
            return tetrahedronIndex;
        }
    }

    // for some reason we were not able to locate the tetrahedron after
    // walking through the structure. now try to find it by looking at
    // every single tetrahedron.
    for(size_t i = 0; i < d->tetrahedra.size(); i++){
        const Tetrahedron &tetrahedron = d->tetrahedra[i];
        if(!tetrahedron.valid){
            continue;
        }

        const Point3 &a = position(tetrahedron.verticies[0]);
        const Point3 &b = position(tetrahedron.verticies[1]);
        const Point3 &c = position(tetrahedron.verticies[2]);
        const Point3 &d = position(tetrahedron.verticies[3]);

        if(chemkit::geometry::planeOrientation(a, b, c, point) < 0 &&
           chemkit::geometry::planeOrientation(a, d, b, point) < 0 &&
           chemkit::geometry::planeOrientation(a, c, d, point) < 0 &&
           chemkit::geometry::planeOrientation(b, d, c, point) < 0){
            return i;
        }
    }

    // should not get here
    return -1;
}

/// Returns a list of tetrahedra that contain the vertex in their
/// circumsphere.
std::vector<int> DelaunayTriangulation::findContainingTetrahedra(int vertex) const
{
    const Point3 &point = position(vertex);

    std::vector<int> tetrahedra;

    int initialTetrahedron = location(point);

    std::set<int> visited;
    std::deque<int> queue;

    queue.push_back(initialTetrahedron);

    while(!queue.empty()){
        int index = queue.front();
        queue.pop_front();
        if(index == -1 || visited.count(index))
            continue;

        visited.insert(index);
        const Tetrahedron &tetrahedron = d->tetrahedra[index];

        int va = tetrahedron.verticies[0];
        int vb = tetrahedron.verticies[1];
        int vc = tetrahedron.verticies[2];
        int vd = tetrahedron.verticies[3];

        Point3 pa = position(va);
        Point3 pb = position(vb);
        Point3 pc = position(vc);
        Point3 pd = position(vd);

        if(chemkit::geometry::planeOrientation(pa, pb, pc, pd) < 0){
            std::swap(pa, pb);
            std::swap(va, vb);
        }

        if(isWeighted()){
            Real wa = weight(va);
            Real wb = weight(vb);
            Real wc = weight(vc);
            Real wd = weight(vd);
            Real wp = weight(vertex);

            if(chemkit::geometry::sphereOrientation(pa, pb, pc, pd, point, wa, wb, wc, wd, wp) > 0){
                tetrahedra.push_back(index);

                for(int i = 0; i < 4; i++){
                    queue.push_back(tetrahedron.neighbors[i]);
                }
            }
        }
        else{
            if(chemkit::geometry::sphereOrientation(pa, pb, pc, pd, point) > 0){
                tetrahedra.push_back(index);

                for(int i = 0; i < 4; i++){
                    queue.push_back(tetrahedron.neighbors[i]);
                }
            }
        }
    }

    return tetrahedra;
}

void DelaunayTriangulation::insertPoint(int index)
{
    Point3 point = position(index);

    const std::vector<int> containingTetrahedra = findContainingTetrahedra(index);

    std::vector<Triangle> faces;
    std::vector<std::pair<int, int> > faceNeighbor;
    std::vector<int> count;

    foreach(int index, containingTetrahedra){
        const Tetrahedron &tetrahedron = d->tetrahedra[index];
        for(int i = 0; i < 4; i++){
            for(int j = i + 1; j < 4; j++){
                for(int k = j + 1; k < 4; k++){
                    Triangle face;
                    face[0] = tetrahedron.verticies[i];
                    face[1] = tetrahedron.verticies[j];
                    face[2] = tetrahedron.verticies[k];
                    std::sort(face.begin(), face.end());

                    unsigned int faceIndex = std::distance(faces.begin(), std::find(faces.begin(), faces.end(), face));
                    if(faceIndex == faces.size()){
                        faces.push_back(face);
                        count.push_back(1);

                        int faceNumber;
                        if(i == 0 && j == 1 && k == 2)
                            faceNumber = 0;
                        else if(i == 0 && j == 1 && k == 3)
                            faceNumber = 1;
                        else if(i == 0 && j == 2 && k == 3)
                            faceNumber = 2;
                        else
                            faceNumber = 3;

                        faceNeighbor.push_back(std::make_pair(index, faceNumber));
                    }
                    else{
                        count[faceIndex]++;
                    }
                }
            }
        }
    }

    // remove containing tetrahedra
    foreach(int tetrahedron, containingTetrahedra){
        d->tetrahedra[tetrahedron].valid = false;
    }

    // add new tetrahedra
    std::vector<int> newTetrahedra;
    for(unsigned int i = 0; i < faces.size(); i++){
        if(count[i] == 1){
            const Triangle &face = faces[i];
            Tetrahedron tetrahedron;
            int tetrahedronIndex = d->tetrahedra.size();

            const Point3 &a = position(face[0]);
            const Point3 &b = position(face[1]);
            const Point3 &c = position(face[2]);

            if(chemkit::geometry::planeOrientation(a, b, c, point) < 0){
                tetrahedron.verticies[0] = face[0];
                tetrahedron.verticies[1] = face[1];
                tetrahedron.verticies[2] = face[2];
                tetrahedron.verticies[3] = index;
            }
            else{
                tetrahedron.verticies[0] = face[0];
                tetrahedron.verticies[1] = face[2];
                tetrahedron.verticies[2] = face[1];
                tetrahedron.verticies[3] = index;
            }

            std::pair<int, int> neighbor = faceNeighbor[i];
            int neighborIndex = neighbor.first;
            int neighborFace = neighbor.second;

            int oldNeighborIndex = d->tetrahedra[neighborIndex].neighbors[neighborFace];
            if(oldNeighborIndex != -1){
                const Tetrahedron &oldNeighbor = d->tetrahedra[oldNeighborIndex];

                int oldNeighborFace = 0;
                for(int j = 0; j < 4; j++){
                    if(oldNeighbor.neighbors[j] == neighborIndex){
                        oldNeighborFace = j;
                    }
                }

                d->tetrahedra[oldNeighborIndex].neighbors[oldNeighborFace] = tetrahedronIndex;
            }

            tetrahedron.neighbors[0] = oldNeighborIndex; // abc
            tetrahedron.neighbors[1] = -2; // abd
            tetrahedron.neighbors[2] = -2; // acd
            tetrahedron.neighbors[3] = -2; // bcd

            tetrahedron.valid = true;
            d->tetrahedra.push_back(tetrahedron);
            newTetrahedra.push_back(tetrahedronIndex);
        }
    }

    // fix up neighbors in new tetrahedra
    for(unsigned int i = 0; i < newTetrahedra.size(); i++){
        Tetrahedron &tetrahedron = d->tetrahedra[newTetrahedra[i]];

        for(unsigned int j = 0; j < newTetrahedra.size(); j++){
            if(i == j){
                continue;
            }

            int otherIndex = newTetrahedra[j];
            const Tetrahedron &other = d->tetrahedra[otherIndex];
            boost::array<int, 4> verticies;
            verticies[0] = other.verticies[0];
            verticies[1] = other.verticies[1];
            verticies[2] = other.verticies[2];
            verticies[3] = other.verticies[3];

            if(std::find(verticies.begin(), verticies.end(), tetrahedron.verticies[0]) != verticies.end() &&
               std::find(verticies.begin(), verticies.end(), tetrahedron.verticies[1]) != verticies.end() &&
               std::find(verticies.begin(), verticies.end(), tetrahedron.verticies[3]) != verticies.end()){
                tetrahedron.neighbors[1] = otherIndex; // abd
            }
            else if(std::find(verticies.begin(), verticies.end(), tetrahedron.verticies[0]) != verticies.end() &&
                    std::find(verticies.begin(), verticies.end(), tetrahedron.verticies[2]) != verticies.end() &&
                    std::find(verticies.begin(), verticies.end(), tetrahedron.verticies[3]) != verticies.end()){
                tetrahedron.neighbors[2] = otherIndex; // acd
            }
            else if(std::find(verticies.begin(), verticies.end(), tetrahedron.verticies[1]) != verticies.end() &&
                    std::find(verticies.begin(), verticies.end(), tetrahedron.verticies[2]) != verticies.end() &&
                    std::find(verticies.begin(), verticies.end(), tetrahedron.verticies[3]) != verticies.end()){
                tetrahedron.neighbors[3] = otherIndex; // bcd
            }
        }
    }
}

bool DelaunayTriangulation::isExternal(int index) const
{
    const Tetrahedron &tetrahedron = d->tetrahedra[index];

    for(int i = 0; i < 4; i++){
        unsigned int vertex = tetrahedron.verticies[i];

        if(vertex >= (d->verticies.size() - 4)){
            return true;
        }
    }

    return false;
}

} // end chemkit namespace
