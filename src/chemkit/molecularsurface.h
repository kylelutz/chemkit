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

#ifndef CHEMKIT_MOLECULARSURFACE_H
#define CHEMKIT_MOLECULARSURFACE_H

#include "chemkit.h"

#include <QFuture>

#include "point.h"

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
        void setProbeRadius(Float radius);
        Float probeRadius() const;

        // geometry
        Point position(int index) const;
        Float radius(int index) const;
        Float volume() const;
        QFuture<Float> volumeAsync() const;
        Float surfaceArea() const;
        QFuture<Float> surfaceAreaAsync() const;

    private:
        // internal methods
        const AlphaShape* alphaShape() const;
        void setCalculated(bool calculated) const;
        Float intersectionArea(int i, int j) const;
        Float intersectionArea(int i, int j, int k) const;
        Float intersectionArea(int i, int j, int k, int l) const;
        Float intersectionVolume(int i, int j) const;
        Float intersectionVolume(int i, int j, int k) const;
        Float intersectionVolume(int i, int j, int k, int l) const;
        Float ballArea(int index) const;
        Float capHeight(int i, int j) const;
        Float capArea(int i, int j) const;
        Float capVolume(int i, int j) const;
        Float cap2Area(int i, int j, int k) const;
        Float cap2Volume(int i, int j, int k) const;
        Float cap3Area(int i, int j, int k, int l) const;
        Float cap3Volume(int i, int j, int k, int l) const;
        Float diskArea(int i, int j) const;
        Float diskLength(int i, int j) const;
        Float diskRadius(int i, int j) const;
        Point triangleDual(int i, int j, int k) const;
        Float segmentArea(int i, int j, int k) const;
        Float segmentAngle(int i, int j, int k) const;
        Float segmentLength(int i, int j, int k) const;
        Float segmentHeight(int i, int j, int k) const;
        Float segment2Area(int i, int j, int k, int l) const;
        Float segment2Angle(int i, int j, int k, int l) const;
        Float segment2Length(int i, int j, int k, int l) const;
        bool ccw(int i, int j, int k, int l) const;

    private:
        MolecularSurfacePrivate* const d;
};

} // end chemkit namespace

#endif // CHEMKIT_MOLECULARSURFACE_H
