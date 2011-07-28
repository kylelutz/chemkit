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

// The formulae for sphere intersection area and volume are
// derived from those presented in: "Measuring Space Filling
// Diagrams and Voids" by Herbert Edelsbrunner and Ping Fu.

#include "molecularsurface.h"

#include <QtConcurrentRun>

#include "foreach.h"
#include "vector3.h"
#include "geometry.h"
#include "molecule.h"
#include "alphashape.h"
#include "delaunaytriangulation.h"

namespace chemkit {

namespace {

const Float pi = chemkit::constants::Pi;

Float angleDihedral(const Point3 &s, const Point3 &t, const Point3 &u, const Point3 &v)
{
    Vector3 mu = (u - s).cross(u - t);
    Vector3 mv = (v - s).cross(v - t);

    Vector3 nu = mu.normalized();
    Vector3 nv = mv.normalized();

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
        std::vector<Point3> points;
        std::vector<Float> radii;
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
/// Provides names for each of the available surface types:
///     - \c VanDerWaals
///     - \c SolventAccessible
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
            d->points.push_back(atom->position());
            d->radii.push_back(atom->vanDerWaalsRadius());
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

const AlphaShape* MolecularSurface::alphaShape() const
{
    if(!d->alphaShape){
        // calculate weights (weight = radius sqaured)
        std::vector<Float> weights(d->points.size());
        for(unsigned int i = 0; i < d->points.size(); i++){
            weights[i] = pow(radius(i), 2);
        }

        d->alphaShape = new AlphaShape(d->points, weights);
    }

    return d->alphaShape;
}

// --- Geometry ------------------------------------------------------------ //
/// Returns the position of the sphere at \p index.
Point3 MolecularSurface::position(int index) const
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
        for(unsigned int i = 0; i < d->points.size(); i++){
            Float r = radius(i);

            d->volume += (4.0/3.0) * pi * r*r*r;
        }

        // subtract volume from each edge
        foreach(const std::vector<int> &edge, alphaShape->edges()){
            d->volume -= intersectionVolume(edge[0], edge[1]);
        }

        // add volume from each triangle
        foreach(const std::vector<int> &triangle, alphaShape->triangles()){
            d->volume += intersectionVolume(triangle[0], triangle[1], triangle[2]);
        }

        // subtract volume from each tetrahedron
        foreach(const std::vector<int> &tetrahedron, alphaShape->tetrahedra()){
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
        for(unsigned int i = 0; i < d->points.size(); i++){
            Float r = radius(i);

            d->surfaceArea += 4.0 * pi * r*r;
        }

        // subtract volume and area from each edge
        foreach(const std::vector<int> &edge, alphaShape->edges()){
            d->surfaceArea -= intersectionArea(edge[0], edge[1]);
        }

        // add volume and area from each triangle
        foreach(const std::vector<int> &triangle, alphaShape->triangles()){
            d->surfaceArea += intersectionArea(triangle[0], triangle[1], triangle[2]);
        }

        // subtract volume and area from each tetrahedron
        foreach(const std::vector<int> &tetrahedron, alphaShape->tetrahedra()){
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
    const Point3 &s = position(i);

    Point3 y = d->alphaShape->orthocenter(i, j);

    // check if vertex i is attached to vertex j
    if(d->alphaShape->vertexAttached(i, j)){
        return radius(i) + chemkit::geometry::distance(s, y);
    }
    else{
        return radius(i) - chemkit::geometry::distance(s, y);
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
    Point3 pjk = triangleDual(i, j, k);

    Float lj = segmentAngle(i, j, k);
    Float lk = segmentAngle(i, k, j);

    const Point3 &s = position(i);
    const Point3 &t = position(j);
    const Point3 &u = position(k);

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
        std::swap(k, l);
    }

    const Point3 &s = position(i);
    const Point3 &t = position(j);
    const Point3 &u = position(k);
    const Point3 &v = position(l);

    Point3 pkj = triangleDual(i, k, j);
    Point3 plk = triangleDual(i, l, k);
    Point3 pjl = triangleDual(i, j, l);

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

Point3 MolecularSurface::triangleDual(int i, int j, int k) const
{
    Point3 y = d->alphaShape->orthocenter(i, j, k);

    const Point3 &s = d->points[i];
    const Point3 &t = d->points[j];
    const Point3 &u = d->points[k];

    Vector3 n = (t - s).cross(u - s);

    Vector3 ys = y - s;

    Float s1 = (ys).dot(n);
    Float s2 = n.dot(n);
    Float s3 = (ys).dot(ys);

    Float r = radius(i);

    Float xi = (-s1 + sqrt(s1*s1 - s3 * s2 + r*r * s2)) / s2;

    return y + (n * xi);
}

Float MolecularSurface::segmentArea(int i, int j, int k) const
{
    Float s = (1.0/2.0) * diskRadius(i, j) * segmentLength(i, j, k);

    Point3 pjk = triangleDual(i, j, k);
    Point3 pkj = triangleDual(i, k, j);

    Float h = diskRadius(i, j) - segmentHeight(i, j, k);
    Float t = (1.0/2.0) * h * chemkit::geometry::distance(pjk, pkj);

    return s - t;
}

Float MolecularSurface::segmentAngle(int i, int j, int k) const
{
    Point3 pjk = triangleDual(i, j, k);

    const Point3 &s = d->points[i];
    const Point3 &t = d->points[j];
    const Point3 &u = d->points[k];

    return 2.0 * angleDihedral(s, t, u, pjk);
}

Float MolecularSurface::segmentLength(int i, int j, int k) const
{
    return segmentAngle(i, j, k) * diskLength(i, j);
}

Float MolecularSurface::segmentHeight(int i, int j, int k) const
{
    Point3 y2 = d->alphaShape->orthocenter(i, j);
    Point3 y3 = d->alphaShape->orthocenter(i, j, k);

    // check if vertex k is attached to the edge (i, j)
    if(d->alphaShape->edgeAttached(i, j, k)){
        return diskRadius(i, j) + chemkit::geometry::distance(y2, y3);
    }
    else{
        return diskRadius(i, j) - chemkit::geometry::distance(y2, y3);
    }
}

Float MolecularSurface::segment2Area(int i, int j, int k, int l) const
{
    if(!ccw(i, j, k, l))
        std::swap(k, l);

    Point3 pjk = triangleDual(i, j, k);
    Point3 pkj = triangleDual(i, k, j);
    Point3 pjl = triangleDual(i, j, l);
    Point3 plj = triangleDual(i, l, j);

    Point3 y = d->alphaShape->orthocenter(i, j, k, l);

    Float hk = segmentHeight(i, j, k);
    Float hl = segmentHeight(i, j, l);

    Float rij = diskRadius(i, j);

    Float s = (1.0/2.0) * rij * segment2Length(i, j, k, l);
    Float tk = (1.0/2.0) * (rij - hk) * chemkit::geometry::distance(pkj, y);
    Float tl = (1.0/2.0) * (rij - hl) * chemkit::geometry::distance(pjl, y);

    return s - tk - tl;
}

Float MolecularSurface::segment2Angle(int i, int j, int k, int l) const
{
    Point3 pjl = triangleDual(i, j, l);
    Point3 pkj = triangleDual(i, k, j);

    const Point3 &s = d->points[i];
    const Point3 &t = d->points[j];
    const Point3 &u = d->points[k];
    const Point3 &v = d->points[l];

    return angleDihedral(s, t, u, pkj) + angleDihedral(s, t, v, pjl) - angleDihedral(s, t, u, v);
}

Float MolecularSurface::segment2Length(int i, int j, int k, int l) const
{
    return segment2Angle(i, j, k, l) * diskLength(i, j);
}

bool MolecularSurface::ccw(int i, int j, int k, int l) const
{
    const Point3 &a = position(i);
    const Point3 &b = position(j);
    const Point3 &c = position(k);
    const Point3 &d = position(l);

    return chemkit::geometry::planeOrientation(a, b, c, d) > 0;
}

} // end chemkit namespace
