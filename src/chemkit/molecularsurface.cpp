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

// The formulae for sphere intersection area and volume are
// derived from those presented in: "Measuring Space Filling
// Diagrams and Voids" by Herbert Edelsbrunner and Ping Fu.

#include "molecularsurface.h"

#include "vector.h"
#include "geometry.h"
#include "molecule.h"
#include "alphashape.h"
#include "staticmatrix.h"
#include "delaunaytriangulation.h"

namespace chemkit {

namespace {

const Float pi = chemkit::constants::Pi;

Float angleDihedral(const Point &s, const Point &t, const Point &u, const Point &v)
{
    Vector mu = (u - s).cross(u - t);
    Vector mv = (v - s).cross(v - t);

    Vector nu = mu.normalized();
    Vector nv = mv.normalized();

    return acos(nu.dot(nv)) / (2.0 * pi);
}

} // end anonymous namespace

// === MolecularSurfacePrivate ============================================= //
class MolecularSurfacePrivate
{
    public:
        const Molecule *molecule;
        MolecularSurface::SurfaceType surfaceType;
        Float probeRadius;
        QVector<Point> points;
        QVector<Float> radii;
        AlphaShape *alphaShape;
        Float volume;
        Float surfaceArea;
        bool volumeCalculated;
        bool surfaceAreaCalculated;
};

// === MolecularSurface ==================================================== //
/// \class MolecularSurface molecularsurface.h chemkit/molecularsurface.h
/// \ingroup chemkit
/// \brief  The MolecularSurface class represents a molecular
///         surface.
///
/// The following example shows how to calculate the solvent
/// accessible surface area of a Protein.
/// \code
/// // create the surface object using the protein molecule
/// chemkit::MolecularSurface surface(protein->molecule());
///
/// // set the surface type to solvent accessible
/// surface.setSurfaceType(chemkit::MolecularSurface::SolventAccessible);
///
/// // set the solvent probe radius to 1.4 angstroms
/// surface.setProbeRadius(1.4);
///
/// // calculate the surface area
/// chemkit::Float area = surface.surfaceArea();
/// \endcode

/// \enum MolecularSurface::SurfaceType
/// Provides names for each of the availble surface types:
///     - \c VanDerWaals
///     - \c SolventAccessible,
///     - \c SolventExcluded

// --- Construction and Destruction ---------------------------------------- //
/// Creates a new molecular surface for \p molecule.
MolecularSurface::MolecularSurface(const Molecule *molecule, SurfaceType type)
    : d(new MolecularSurfacePrivate)
{
    d->molecule = molecule;
    d->surfaceType = type;
    d->probeRadius = 1.4;

    if(molecule){
        foreach(const Atom *atom, molecule->atoms()){
            d->points.append(atom->position());
            d->radii.append(atom->vanDerWaalsRadius());
        }
    }

    d->alphaShape = 0;
    d->volumeCalculated = false;
    d->surfaceAreaCalculated = false;
}

/// Destroys the molecular surface object.
MolecularSurface::~MolecularSurface()
{
    delete d->alphaShape;
    delete d;
}

// --- Properties ---------------------------------------------------------- //
/// Sets the molecule for the surface.
void MolecularSurface::setMolecule(const Molecule *molecule)
{
    d->molecule = molecule;

    // update atom positions and radii
    if(molecule){
        d->points.resize(molecule->size());
        d->radii.resize(molecule->size());

        for(int i = 0; i < molecule->size(); i++){
            const Atom *atom = molecule->atom(i);

            d->points[i] = atom->position();
            d->radii[i] = atom->vanDerWaalsRadius();
        }
    }

    setCalculated(false);
}

/// Returns the molecule for the surface.
const Molecule* MolecularSurface::molecule() const
{
    return d->molecule;
}

/// Sets the surface type to \p type.
void MolecularSurface::setSurfaceType(SurfaceType type)
{
    d->surfaceType = type;

    setCalculated(false);
}

/// Returns the surface type.
MolecularSurface::SurfaceType MolecularSurface::surfaceType() const
{
    return d->surfaceType;
}

/// Sets the probe radius to \p radius.
void MolecularSurface::setProbeRadius(Float radius)
{
    d->probeRadius = radius;

    setCalculated(false);
}

/// Returns the probe radius.
///
/// The default probe radius is 1.4 Angstroms which approximates
/// the radius of a water molecule.
Float MolecularSurface::probeRadius() const
{
    return d->probeRadius;
}

// --- Geometry ------------------------------------------------------------ //
/// Returns the position of the sphere at \p index.
Point MolecularSurface::position(int index) const
{
    return d->points[index];
}

/// Returns the radius of the sphere at \p index.
Float MolecularSurface::radius(int index) const
{
    if(d->surfaceType == VanDerWaals)
        return d->radii[index];
    else
        return d->radii[index] + d->probeRadius;
}

/// Returns the total volume of the surface. The returned volume
/// is in Angstroms cubed (\f$ \AA^{3} \f$).
Float MolecularSurface::volume() const
{
    if(!d->volumeCalculated){
        d->volume = 0;

        const AlphaShape *alphaShape = this->alphaShape();

        // add volume and area for each vertex
        for(int i = 0; i < d->points.size(); i++){
            Float r = radius(i);

            d->volume += (4.0/3.0) * pi * r*r*r;
        }

        // subtract volume from each edge
        foreach(QVector<int> edge, alphaShape->edges()){
            d->volume -= intersectionVolume(edge[0], edge[1]);
        }

        // add volume from each triangle
        foreach(const QVector<int> &triangle, alphaShape->triangles()){
            d->volume += intersectionVolume(triangle[0], triangle[1], triangle[2]);
        }

        // subtract volume from each tetrahedron
        foreach(QVector<int> tetrahedron, alphaShape->tetrahedra()){
            d->volume -= intersectionVolume(tetrahedron[0], tetrahedron[1], tetrahedron[2], tetrahedron[3]);
        }

        d->volumeCalculated = true;
    }

    return d->volume;
}

QFuture<Float> MolecularSurface::volumeAsync() const
{
    QFuture<Float> future = QtConcurrent::run(this, &MolecularSurface::volume);

    return future;
}

/// Returns the total surface area of the surface. The returned
/// area is in Angstroms squared (\f$ \AA^{2} \f$).
Float MolecularSurface::surfaceArea() const
{
    if(!d->surfaceAreaCalculated){
        d->surfaceArea = 0;

        const AlphaShape *alphaShape = this->alphaShape();

        // add volume and area for each vertex
        for(int i = 0; i < d->points.size(); i++){
            Float r = radius(i);

            d->surfaceArea += 4.0 * pi * r*r;
        }

        // subtract volume and area from each edge
        foreach(QVector<int> edge, alphaShape->edges()){
            d->surfaceArea -= intersectionArea(edge[0], edge[1]);
        }

        // add volume and area from each triangle
        foreach(const QVector<int> &triangle, alphaShape->triangles()){
            d->surfaceArea += intersectionArea(triangle[0], triangle[1], triangle[2]);
        }

        // subtract volume and area from each tetrahedron
        foreach(QVector<int> tetrahedron, alphaShape->tetrahedra()){
            d->surfaceArea -= intersectionArea(tetrahedron[0], tetrahedron[1], tetrahedron[2], tetrahedron[3]);
        }

        d->surfaceAreaCalculated = true;
    }

    return d->surfaceArea;
}

QFuture<Float> MolecularSurface::surfaceAreaAsync() const
{
    QFuture<Float> future = QtConcurrent::run(this, &MolecularSurface::surfaceArea);

    return future;
}

// --- Internal Methods ---------------------------------------------------- //
const AlphaShape* MolecularSurface::alphaShape() const
{
    if(!d->alphaShape){
        // calculate weights (weight = radius sqaured)
        QVector<Float> weights(d->points.size());
        for(int i = 0; i < d->points.size(); i++){
            weights[i] = pow(radius(i), 2);
        }

        d->alphaShape = new AlphaShape(d->points, weights);
    }

    return d->alphaShape;
}

void MolecularSurface::setCalculated(bool calculated) const
{
    if(calculated == false){
        delete d->alphaShape;
        d->alphaShape = 0;
        d->volumeCalculated = false;
        d->surfaceAreaCalculated = false;
    }
}

/// Returns the area of intersection between spheres \p i and \p j.
Float MolecularSurface::intersectionArea(int i, int j) const
{
    return 2.0 * pi * (radius(i) * capHeight(i, j) + radius(j) * capHeight(j, i));
}

/// Returns the area of intersection between spheres \p i, \p j and
/// \p k.
Float MolecularSurface::intersectionArea(int i, int j, int k) const
{
    return cap2Area(i, j, k) + cap2Area(j, i, k) + cap2Area(k, i, j);
}

/// Returns the area of intersection between spheres \p i, \p j, \p k
/// and \p l.
Float MolecularSurface::intersectionArea(int i, int j, int k, int l) const
{
    return cap3Area(i, j, k, l) + cap3Area(j, i, k, l) + cap3Area(k, i, j, l) + cap3Area(l, i, j, k);
}

/// Returns the volume of intersection between spheres \p i, and
/// \p j.
Float MolecularSurface::intersectionVolume(int i, int j) const
{
    return capVolume(i, j) + capVolume(j, i);
}

/// Returns the volume of intersection between spheres \p i, \p j,
/// and \p k.
Float MolecularSurface::intersectionVolume(int i, int j, int k) const
{
    return cap2Volume(i, j, k) + cap2Volume(j, i, k) + cap2Volume(k, i, j);
}

/// Returns the volume of intersection between spheres \p i, \p j,
/// \p k and \p l.
Float MolecularSurface::intersectionVolume(int i, int j, int k, int l) const
{
    return cap3Volume(i, j, k, l) + cap3Volume(j, i, k, l) + cap3Volume(k, i, j, l) + cap3Volume(l, i, j, k);
}

Float MolecularSurface::ballArea(int index) const
{
    Float r = radius(index);

    return 4.0 * pi * r*r;
}

Float MolecularSurface::capHeight(int i, int j) const
{
    const Point &s = position(i);

    Point y = d->alphaShape->orthocenter(i, j);

    // check if vertex i is attached to vertex j
    if(d->alphaShape->vertexAttached(i, j)){
        return radius(i) + s.distance(y);
    }
    else{
        return radius(i) - s.distance(y);
    }
}

Float MolecularSurface::capArea(int i, int j) const
{
    return 2.0 * pi * radius(i) * capHeight(i, j);
}

Float MolecularSurface::capVolume(int i, int j) const
{
    Float s = radius(i) * capArea(i, j);
    Float c = (radius(i) - capHeight(i, j)) * diskArea(i, j);

    return (1.0/3.0) * (s - c);
}

Float MolecularSurface::cap2Area(int i, int j, int k) const
{
    Point pjk = triangleDual(i, j, k);

    Float lj = segmentAngle(i, j, k);
    Float lk = segmentAngle(i, k, j);

    const Point &s = position(i);
    const Point &t = position(j);
    const Point &u = position(k);

    Float r = radius(i);
    Float phi = (1.0/2.0) - angleDihedral(s, pjk, t, u);

    Float a1 = ballArea(i) * phi;
    Float a2 = 2.0 * pi * r * lj * (r - capHeight(i, j));
    Float a3 = 2.0 * pi * r * lk * (r - capHeight(i, k));

    return a1 - a2 - a3;
}

Float MolecularSurface::cap2Volume(int i, int j, int k) const
{
    Float s2 = (1.0/3.0) * radius(i) * cap2Area(i, j, k);
    Float cj = (1.0/3.0) * (radius(i) - capHeight(i, j)) * segmentArea(i, j, k);
    Float ck = (1.0/3.0) * (radius(i) - capHeight(i, k)) * segmentArea(i, k, j);

    return s2 - cj - ck;
}

Float MolecularSurface::cap3Area(int i, int j, int k, int l) const
{
    if(!ccw(i, j, k, l)){
        qSwap(k, l);
    }

    const Point &s = position(i);
    const Point &t = position(j);
    const Point &u = position(k);
    const Point &v = position(l);

    Point pkj = triangleDual(i, k, j);
    Point plk = triangleDual(i, l, k);
    Point pjl = triangleDual(i, j, l);

    Float lj = segment2Angle(i, j, k, l);
    Float lk = segment2Angle(i, k, l, j);
    Float ll = segment2Angle(i, l, j, k);

    Float rho_kj = (1.0/2.0) - angleDihedral(s, pkj, u, t);
    Float rho_lk = (1.0/2.0) - angleDihedral(s, plk, v, u);
    Float rho_jl = (1.0/2.0) - angleDihedral(s, pjl, t, v);

    Float a1 = (1.0/2.0) * ballArea(i) * (rho_kj + rho_lk + rho_jl - (1.0/2.0));
    Float a2 = 2.0 * pi * radius(i) * lj * (radius(i) - capHeight(i, j));
    Float a3 = 2.0 * pi * radius(i) * lk * (radius(i) - capHeight(i, k));
    Float a4 = 2.0 * pi * radius(i) * ll * (radius(i) - capHeight(i, l));

    return a1 - a2 - a3 - a4;
}

Float MolecularSurface::cap3Volume(int i, int j, int k, int l) const
{
    Float s3 = (1.0/3.0) * radius(i) * cap3Area(i, j, k, l);
    Float cj = (1.0/3.0) * (radius(i) - capHeight(i, j)) * segment2Area(i, j, k, l);
    Float ck = (1.0/3.0) * (radius(i) - capHeight(i, k)) * segment2Area(i, k, j, l);
    Float cl = (1.0/3.0) * (radius(i) - capHeight(i, l)) * segment2Area(i, l, j, k);

    return s3 - cj - ck - cl;
}

Float MolecularSurface::diskArea(int i, int j) const
{
    return (1.0/2.0) * diskRadius(i, j) * diskLength(i, j);
}

Float MolecularSurface::diskLength(int i, int j) const
{
    return 2.0 * pi * diskRadius(i, j);
}

Float MolecularSurface::diskRadius(int i, int j) const
{
    return sqrt(capHeight(i, j) * (2.0 * radius(i) - capHeight(i, j)));
}

Point MolecularSurface::triangleDual(int i, int j, int k) const
{
    Point y = d->alphaShape->orthocenter(i, j, k);

    const Point &s = d->points[i];
    const Point &t = d->points[j];
    const Point &u = d->points[k];

    Vector n = (t - s).cross(u - s);

    Vector ys = y - s;

    Float s1 = (ys).dot(n);
    Float s2 = n.dot(n);
    Float s3 = (ys).dot(ys);

    Float r = radius(i);

    Float xi = (-s1 + sqrt(s1*s1 - s3 * s2 + r*r * s2)) / s2;

    return y + (n.scaled(xi));
}

Float MolecularSurface::segmentArea(int i, int j, int k) const
{
    Float s = (1.0/2.0) * diskRadius(i, j) * segmentLength(i, j, k);

    Point pjk = triangleDual(i, j, k);
    Point pkj = triangleDual(i, k, j);

    Float h = diskRadius(i, j) - segmentHeight(i, j, k);
    Float t = (1.0/2.0) * h * pjk.distance(pkj);

    return s - t;
}

Float MolecularSurface::segmentAngle(int i, int j, int k) const
{
    Point pjk = triangleDual(i, j, k);

    const Point &s = d->points[i];
    const Point &t = d->points[j];
    const Point &u = d->points[k];

    return 2.0 * angleDihedral(s, t, u, pjk);
}

Float MolecularSurface::segmentLength(int i, int j, int k) const
{
    return segmentAngle(i, j, k) * diskLength(i, j);
}

Float MolecularSurface::segmentHeight(int i, int j, int k) const
{
    Point y2 = d->alphaShape->orthocenter(i, j);
    Point y3 = d->alphaShape->orthocenter(i, j, k);

    // check if vertex k is attached to the edge (i, j)
    if(d->alphaShape->edgeAttached(i, j, k)){
        return diskRadius(i, j) + y2.distance(y3);
    }
    else{
        return diskRadius(i, j) - y2.distance(y3);
    }
}

Float MolecularSurface::segment2Area(int i, int j, int k, int l) const
{
    if(!ccw(i, j, k, l))
        qSwap(k, l);

    Point pjk = triangleDual(i, j, k);
    Point pkj = triangleDual(i, k, j);
    Point pjl = triangleDual(i, j, l);
    Point plj = triangleDual(i, l, j);

    Point y = d->alphaShape->orthocenter(i, j, k, l);

    Float hk = segmentHeight(i, j, k);
    Float hl = segmentHeight(i, j, l);

    Float rij = diskRadius(i, j);

    Float s = (1.0/2.0) * rij * segment2Length(i, j, k, l);
    Float tk = (1.0/2.0) * (rij - hk) * pkj.distance(y);
    Float tl = (1.0/2.0) * (rij - hl) * pjl.distance(y);

    return s - tk - tl;
}

Float MolecularSurface::segment2Angle(int i, int j, int k, int l) const
{
    Point pjl = triangleDual(i, j, l);
    Point pkj = triangleDual(i, k, j);

    const Point &s = d->points[i];
    const Point &t = d->points[j];
    const Point &u = d->points[k];
    const Point &v = d->points[l];

    return angleDihedral(s, t, u, pkj) + angleDihedral(s, t, v, pjl) - angleDihedral(s, t, u, v);
}

Float MolecularSurface::segment2Length(int i, int j, int k, int l) const
{
    return segment2Angle(i, j, k, l) * diskLength(i, j);
}

bool MolecularSurface::ccw(int i, int j, int k, int l) const
{
    const Point &a = position(i);
    const Point &b = position(j);
    const Point &c = position(k);
    const Point &d = position(l);

    return chemkit::geometry::planeOrientation(a, b, c, d) > 0;
}

} // end chemkit namespace
