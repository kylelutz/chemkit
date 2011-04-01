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

#ifndef CHEMKIT_ALPHASHAPE_H
#define CHEMKIT_ALPHASHAPE_H

#include <vector>

#include <QtCore>

#include "point3.h"

namespace chemkit {

class AlphaShapePrivate;

class CHEMKIT_EXPORT AlphaShape
{
    public:
        // enumerations
        enum Classification {
            Interior = 0x01,
            Regular = 0x02,
            Singular = 0x04
        };

        // construction and destruction
        AlphaShape(const std::vector<Point3> &points);
        AlphaShape(const std::vector<Point3> &points, const std::vector<Float> &weights);
        ~AlphaShape();

        // properties
        int size() const;
        Point3 position(int vertex) const;
        Float weight(int vertex) const;
        void setAlphaValue(Float alphaValue);
        Float alphaValue() const;

        // simplicies
        QList<int> verticies() const;
        int vertexCount() const;
        QList<std::vector<int> > edges() const;
        int edgeCount() const;
        QList<std::vector<int> > triangles() const;
        int triangleCount() const;
        QList<std::vector<int> > tetrahedra() const;
        int tetrahedronCount() const;

        // geometry
        Float volume() const;
        Float surfaceArea() const;
        Point3 orthocenter(int i, int j) const;
        Point3 orthocenter(int i, int j, int k) const;
        Point3 orthocenter(int i, int j, int k, int l) const;
        Float orthoradius(int i, int j) const;
        Float orthoradius(int i, int j, int k) const;
        Float orthoradius(int i, int j, int k, int l) const;
        bool vertexAttached(int i, int j) const;
        bool edgeAttached(int i, int j, int k) const;
        bool triangleAttached(int i, int j, int k, int l) const;
        bool triangleAttached(int i, int j, int k, int l, int m) const;

    private:
        AlphaShapePrivate* const d;
};

} // end chemkit namespace

#endif // CHEMKIT_ALPHASHAPE_H
