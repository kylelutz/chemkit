/******************************************************************************
**
** Copyright (C) 2009-2011 Kyle Lutz <kyle.r.lutz@gmail.com>
**
** This file is part of chemkit. For more information see
** <http://www.chemkit.org>.
**
** chemkit is free software: you can redistribute it and/or modify
** it under the terms of the GNU Lesser General Public License as published by
** the Free Software Foundation, either version 3 of the License, or
** (at your option) any later version.
**
** chemkit is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU Lesser General Public License for more details.
**
** You should have received a copy of the GNU Lesser General Public License
** along with chemkit. If not, see <http://www.gnu.org/licenses/>.
**
******************************************************************************/

#include "delaunaytriangulation.h"

#include "point3.h"
#include "vector3.h"
#include "geometry.h"
#include "alphashape.h"
#include "staticmatrix.h"

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
        QVector<QSet<int> > m_edges;
};

EdgeSet::EdgeSet(int vertexCount)
    : m_edges(vertexCount)
{
}

void EdgeSet::insert(int a, int b)
{
    if(a > b)
        qSwap(a, b);

    m_edges[a].insert(b);
}

bool EdgeSet::contains(int a, int b)
{
    if(a > b)
        qSwap(a, b);

    return m_edges[a].contains(b);
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
        QVector<int> triangle(int index) const;
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

QVector<int> Tetrahedron::triangle(int index) const
{
    QVector<int> triangle(3);

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
        QVector<Point3> verticies;
        QVector<Float> weights;
        QVector<Tetrahedron> tetrahedra;

        bool alphaShapeCalculated;

        QList<QVector<int> > delaunayEdges;
        QList<QVector<int> > delaunayTriangles;
        QList<QVector<int> > delaunayTetrahedra;

        QList<QVector<int> > alphaShapeEdges;
        QList<QVector<int> > alphaShapeTriangles;
        QList<QVector<int> > alphaShapeTetrahedra;
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
DelaunayTriangulation::DelaunayTriangulation(const QVector<Point3> &points)
    : d(new DelaunayTriangulationPrivate)
{
    d->verticies = points;

    d->alphaShapeCalculated = false;

    triangulate(false);
}

/// Creates a new weighted delaunay triangulation for \p points with
/// \p weights.
DelaunayTriangulation::DelaunayTriangulation(const QVector<Point3> &points, const QVector<Float> &weights)
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
Float DelaunayTriangulation::weight(int vertex) const
{
    return d->weights[vertex];
}

/// Returns \c true if the delaunay triangulation has weighted
/// verticies.
bool DelaunayTriangulation::isWeighted() const
{
    return !d->weights.isEmpty();
}

// --- Simplicies ---------------------------------------------------------- //
/// Returns a list of verticies in the delaunay triangulation.
QList<int> DelaunayTriangulation::verticies() const
{
    QList<int> verticies;

    for(int i = 0; i < d->verticies.size() - 4; i++){
        verticies.append(i);
    }

    return verticies;
}

/// Returns the number of verticies in the delaunay triangulation.
int DelaunayTriangulation::vertexCount() const
{
    return verticies().size();
}

/// Returns a list of edges in the delaunay triangulation.
QList<QVector<int> > DelaunayTriangulation::edges() const
{
    if(d->delaunayEdges.isEmpty()){
        QList<QVector<int> > edges;
        EdgeSet edgeSet(d->verticies.size());

        foreach(const QVector<int> &tetrahedron, tetrahedra()){
            for(int i = 0; i < 4; i++){
                for(int j = i + 1; j < 4; j++){
                    QVector<int> edge(2);
                    edge[0] = tetrahedron[i];
                    edge[1] = tetrahedron[j];

                    if(!edgeSet.contains(edge[0], edge[1])){
                        edges.append(edge);
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
QList<QVector<int> > DelaunayTriangulation::triangles() const
{
    if(d->delaunayTriangles.isEmpty()){
        int initialTetrahedron = 0;
        for(int i = 0; i < d->tetrahedra.size(); i++){
            const Tetrahedron &tetrahedron = d->tetrahedra[i];
            if(!tetrahedron.valid || isExternal(i)){
                continue;
            }

            initialTetrahedron = i;
            break;
        }

        QSet<int> visited;
        QStack<int> stack;

        stack.push(initialTetrahedron);

        while(!stack.isEmpty()){
            int index = stack.pop();
            visited.insert(index);
            const Tetrahedron &tetrahedron = d->tetrahedra[index];

            for(int i = 0; i < 4; i++){
                int neighborIndex = tetrahedron.neighbors[i];
                if(neighborIndex == -1 || visited.contains(neighborIndex)){
                    continue;
                }

                d->delaunayTriangles.append(tetrahedron.triangle(i));

                stack.push(neighborIndex);
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
QList<QVector<int> > DelaunayTriangulation::tetrahedra() const
{
    if(d->delaunayTetrahedra.isEmpty()){
        QList<QVector<int> > tetrahedra;

        foreach(const Tetrahedron &tetrahedron, d->tetrahedra){
            if(!tetrahedron.valid){
                continue;
            }

            QVector<int> verticies(4);
            bool external = false;
            for(int i = 0; i < 4; i++){
                int vertex = tetrahedron.verticies[i];

                if(vertex >= (d->verticies.size() - 4)){
                    external = true;
                    break;
                }

                verticies[i] = vertex;
            }

            if(external){
                continue;
            }

            tetrahedra.append(verticies);
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
Float DelaunayTriangulation::volume() const
{
    Float volume = 0;

    foreach(const QVector<int> &tetrahedron, tetrahedra()){
        const Point3 &a = position(tetrahedron[0]);
        const Point3 &b = position(tetrahedron[1]);
        const Point3 &c = position(tetrahedron[2]);
        const Point3 &d = position(tetrahedron[3]);

        volume += chemkit::geometry::tetrahedronVolume(a, b, c, d);
    }

    return volume;
}

/// Returns the total surface area of the triangulation.
Float DelaunayTriangulation::surfaceArea() const
{
    return 0;
}

// --- Alpha Shape --------------------------------------------------------- //
QList<QVector<int> > DelaunayTriangulation::alphaShapeEdges(const AlphaShape *alphaShape) const
{
    if(d->alphaShapeEdges.isEmpty()){
        QList<QVector<int> > alphaShapeEdges;
        EdgeSet alphaEdgeSet(d->verticies.size());

        foreach(const QVector<int> &triangle, alphaShapeTriangles(alphaShape)){
            for(int i = 0; i < 3; i++){
                for(int j = i + 1; j < 3; j++){
                    QVector<int> edge(2);
                    edge[0] = triangle[i];
                    edge[1] = triangle[j];

                    if(!alphaEdgeSet.contains(edge[0], edge[1])){
                        alphaShapeEdges.append(edge);
                        alphaEdgeSet.insert(edge[0], edge[1]);
                    }

                }
            }
        }

        EdgeSet attachedEdgeSet(d->verticies.size());

        foreach(const QVector<int> &triangle, triangles()){
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

        foreach(const QVector<int> edge, edges()){
            if(alphaEdgeSet.contains(edge[0], edge[1])){
                continue;
            }
            else if(attachedEdgeSet.contains(edge[0], edge[1])){
                continue;
            }

            if(alphaShape->orthoradius(edge[0], edge[1]) < alphaShape->alphaValue()){
                alphaShapeEdges.append(edge);
            }
        }

        d->alphaShapeEdges = alphaShapeEdges;
    }

    return d->alphaShapeEdges;
}

QList<QVector<int> > DelaunayTriangulation::alphaShapeTriangles(const AlphaShape *alphaShape) const
{
    if(d->alphaShapeTriangles.isEmpty()){
        calculateAlphaShape(alphaShape);

        QList<QVector<int> > triangles;

        int initialTetrahedron = 0;
        for(int i = 0; i < d->tetrahedra.size(); i++){
            const Tetrahedron &tetrahedron = d->tetrahedra[i];
            if(!tetrahedron.valid || isExternal(i)){
                continue;
            }

            initialTetrahedron = i;
            break;
        }

        QSet<int> visited;
        QStack<int> stack;

        stack.push(initialTetrahedron);

        while(!stack.isEmpty()){
            int index = stack.pop();
            visited.insert(index);
            const Tetrahedron &tetrahedron = d->tetrahedra[index];

            for(int triangleIndex = 0; triangleIndex < 4; triangleIndex++){
                int neighborIndex = tetrahedron.neighbors[triangleIndex];
                if(neighborIndex == -1 || visited.contains(neighborIndex)){
                    continue;
                }

                const Tetrahedron &neighbor = d->tetrahedra[neighborIndex];
                bool neighborExternal = isExternal(neighborIndex);

                if(!neighborExternal){
                    stack.push(neighborIndex);
                }
                else{
                    visited.insert(neighborIndex);
                }

                QVector<int> triangle = tetrahedron.triangle(triangleIndex);

                if(tetrahedron.inAlphaShape && (!neighborExternal && neighbor.inAlphaShape)){
                    triangles.append(triangle);
                }
                else if(tetrahedron.inAlphaShape || (!neighborExternal && neighbor.inAlphaShape)){
                    triangles.append(triangle);
                }
                else{
                    QList<int> tetrahedronVerticies;
                    QList<int> neighborVerticies;
                    for(int i = 0; i < 4; i++){
                        tetrahedronVerticies.append(tetrahedron.verticies[i]);
                        neighborVerticies.append(neighbor.verticies[i]);
                    }

                    for(int i = 0; i < 3; i++){
                        tetrahedronVerticies.removeOne(triangle[i]);
                        neighborVerticies.removeOne(triangle[i]);
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
                        triangles.append(triangle);
                    }
                }
            }
        }

        d->alphaShapeTriangles = triangles;
    }

    return d->alphaShapeTriangles;
}

QList<QVector<int> > DelaunayTriangulation::alphaShapeTetrahedra(const AlphaShape *alphaShape) const
{
    if(d->alphaShapeTetrahedra.isEmpty()){
        calculateAlphaShape(alphaShape);

        for(int i = 0; i < d->tetrahedra.size(); i++){
            const Tetrahedron &tetrahedron = d->tetrahedra[i];
            if(!tetrahedron.valid || isExternal(i))
                continue;

            if(tetrahedron.inAlphaShape){
                QVector<int> verticies(4);
                for(int j = 0; j < 4; j++){
                    verticies[j] = tetrahedron.verticies[j];
                }

                d->alphaShapeTetrahedra.append(verticies);
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

    for(int i = 0; i < d->tetrahedra.size(); i++){
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
    d->verticies.append(Point3(0, 1e10, 0));
    d->verticies.append(Point3(1e10, -1e10, 1e10));
    d->verticies.append(Point3(-1e10, -1e10, 1e10));
    d->verticies.append(Point3(0, -1e10, -1e10));

    if(weighted){
        d->weights.append(0);
        d->weights.append(0);
        d->weights.append(0);
        d->weights.append(0);
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
    d->tetrahedra.append(big);

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

    for(;;){
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
            break;
        }
    }

    return tetrahedronIndex;
}

/// Returns a list of tetrahedra that contain the vertex in their
/// circumsphere.
QList<int> DelaunayTriangulation::findContainingTetrahedra(int vertex) const
{
    const Point3 &point = position(vertex);

    QList<int> tetrahedra;

    int initialTetrahedron = location(point);

    QSet<int> visited;
    QQueue<int> queue;

    queue.enqueue(initialTetrahedron);

    while(!queue.isEmpty()){
        int index = queue.dequeue();
        if(index == -1 || visited.contains(index))
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
            qSwap(pa, pb);
            qSwap(va, vb);
        }

        if(isWeighted()){
            Float wa = weight(va);
            Float wb = weight(vb);
            Float wc = weight(vc);
            Float wd = weight(vd);
            Float wp = weight(vertex);

            if(chemkit::geometry::sphereOrientation(pa, pb, pc, pd, point, wa, wb, wc, wd, wp) > 0){
                tetrahedra.append(index);

                for(int i = 0; i < 4; i++){
                    queue.enqueue(tetrahedron.neighbors[i]);
                }
            }
        }
        else{
            if(chemkit::geometry::sphereOrientation(pa, pb, pc, pd, point) > 0){
                tetrahedra.append(index);

                for(int i = 0; i < 4; i++){
                    queue.enqueue(tetrahedron.neighbors[i]);
                }
            }
        }
    }

    return tetrahedra;
}

void DelaunayTriangulation::insertPoint(int index)
{
    Point3 point = position(index);

    QList<int> containingTetrahedra = findContainingTetrahedra(index);

    QList<QVector<int> > faces;
    QList<QPair<int, int> > faceNeighbor;
    QList<int> count;

    foreach(int index, containingTetrahedra){
        const Tetrahedron &tetrahedron = d->tetrahedra[index];
        for(int i = 0; i < 4; i++){
            for(int j = i + 1; j < 4; j++){
                for(int k = j + 1; k < 4; k++){
                    QVector<int> face(3);
                    face[0] = tetrahedron.verticies[i];
                    face[1] = tetrahedron.verticies[j];
                    face[2] = tetrahedron.verticies[k];
                    qSort(face);

                    int faceIndex = faces.indexOf(face);
                    if(faceIndex == -1){
                        faces.append(face);
                        count.append(1);

                        int faceNumber;
                        if(i == 0 && j == 1 && k == 2)
                            faceNumber = 0;
                        else if(i == 0 && j == 1 && k == 3)
                            faceNumber = 1;
                        else if(i == 0 && j == 2 && k == 3)
                            faceNumber = 2;
                        else
                            faceNumber = 3;

                        faceNeighbor.append(qMakePair(index, faceNumber));
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
    QList<int> newTetrahedra;
    for(int i = 0; i < faces.size(); i++){
        if(count[i] == 1){
            QVector<int> face = faces[i];
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

            QPair<int, int> neighbor = faceNeighbor[i];
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
            d->tetrahedra.append(tetrahedron);
            newTetrahedra.append(tetrahedronIndex);
        }
    }

    // fix up neighbors in new tetrahedra
    for(int i = 0; i < newTetrahedra.size(); i++){
        Tetrahedron &tetrahedron = d->tetrahedra[newTetrahedra[i]];

        for(int j = 0; j < newTetrahedra.size(); j++){
            if(i == j){
                continue;
            }

            int otherIndex = newTetrahedra[j];
            const Tetrahedron &other = d->tetrahedra[otherIndex];
            QVector<int> verticies(4);
            verticies[0] = other.verticies[0];
            verticies[1] = other.verticies[1];
            verticies[2] = other.verticies[2];
            verticies[3] = other.verticies[3];

            if(verticies.contains(tetrahedron.verticies[0]) &&
               verticies.contains(tetrahedron.verticies[1]) &&
               verticies.contains(tetrahedron.verticies[3])){
                tetrahedron.neighbors[1] = otherIndex; // abd
            }
            else if(verticies.contains(tetrahedron.verticies[0]) &&
                    verticies.contains(tetrahedron.verticies[2]) &&
                    verticies.contains(tetrahedron.verticies[3])){
                tetrahedron.neighbors[2] = otherIndex; // acd
            }
            else if(verticies.contains(tetrahedron.verticies[1]) &&
                    verticies.contains(tetrahedron.verticies[2]) &&
                    verticies.contains(tetrahedron.verticies[3])){
                tetrahedron.neighbors[3] = otherIndex; // bcd
            }
        }
    }
}

bool DelaunayTriangulation::isExternal(int index) const
{
    const Tetrahedron &tetrahedron = d->tetrahedra[index];

    for(int i = 0; i < 4; i++){
        int vertex = tetrahedron.verticies[i];

        if(vertex >= (d->verticies.size() - 4)){
            return true;
        }
    }

    return false;
}

} // end chemkit namespace
