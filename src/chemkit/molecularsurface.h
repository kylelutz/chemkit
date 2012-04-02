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

#ifndef CHEMKIT_MOLECULARSURFACE_H
#define CHEMKIT_MOLECULARSURFACE_H

#include "chemkit.h"

#include <boost/thread/future.hpp>

#include "point3.h"

namespace chemkit {

class Molecule;
class AlphaShape;
class MolecularSurfacePrivate;

class CHEMKIT_EXPORT MolecularSurface
{
public:
    // enumerations
    enum SurfaceType {
        VanDerWaals,
        SolventAccessible,
        SolventExcluded
    };

    // construction and destruction
    MolecularSurface(const Molecule *molecule = 0, SurfaceType type = VanDerWaals);
    ~MolecularSurface();

    // properties
    void setMolecule(const Molecule *molecule);
    const Molecule* molecule() const;
    void setSurfaceType(SurfaceType type);
    SurfaceType surfaceType() const;
    void setProbeRadius(Real radius);
    Real probeRadius() const;
    const AlphaShape* alphaShape() const;

    // geometry
    Point3 position(int index) const;
    Real radius(int index) const;
    Real volume() const;
    boost::shared_future<Real> volumeAsync() const;
    Real surfaceArea() const;
    boost::shared_future<Real> surfaceAreaAsync() const;

private:
    // internal methods
    void setCalculated(bool calculated) const;
    Real intersectionArea(int i, int j) const;
    Real intersectionArea(int i, int j, int k) const;
    Real intersectionArea(int i, int j, int k, int l) const;
    Real intersectionVolume(int i, int j) const;
    Real intersectionVolume(int i, int j, int k) const;
    Real intersectionVolume(int i, int j, int k, int l) const;
    Real ballArea(int index) const;
    Real capHeight(int i, int j) const;
    Real capArea(int i, int j) const;
    Real capVolume(int i, int j) const;
    Real cap2Area(int i, int j, int k) const;
    Real cap2Volume(int i, int j, int k) const;
    Real cap3Area(int i, int j, int k, int l) const;
    Real cap3Volume(int i, int j, int k, int l) const;
    Real diskArea(int i, int j) const;
    Real diskLength(int i, int j) const;
    Real diskRadius(int i, int j) const;
    Point3 triangleDual(int i, int j, int k) const;
    Real segmentArea(int i, int j, int k) const;
    Real segmentAngle(int i, int j, int k) const;
    Real segmentLength(int i, int j, int k) const;
    Real segmentHeight(int i, int j, int k) const;
    Real segment2Area(int i, int j, int k, int l) const;
    Real segment2Angle(int i, int j, int k, int l) const;
    Real segment2Length(int i, int j, int k, int l) const;
    bool ccw(int i, int j, int k, int l) const;

private:
    MolecularSurfacePrivate* const d;
};

} // end chemkit namespace

#endif // CHEMKIT_MOLECULARSURFACE_H
