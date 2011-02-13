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

#ifndef CHEMKIT_GRAPHICSVECTOR_H
#define CHEMKIT_GRAPHICSVECTOR_H

#include "graphics.h"

#include <chemkit/vector.h>
#include <chemkit/genericvector.h>

namespace chemkit {

class CHEMKIT_GRAPHICS_EXPORT GraphicsVector : public GenericVector<GraphicsFloat>
{
    public:
        // construction and destruction
        GraphicsVector();
        GraphicsVector(GraphicsFloat x, GraphicsFloat y, GraphicsFloat z);
        GraphicsVector(const GenericVector<float> &vector);
        GraphicsVector(const GenericVector<double> &vector);
        GraphicsVector(const StaticVector<float, 3> &vector);
        GraphicsVector(const StaticVector<double, 3> &vector);

        // properties
        Vector toVector() const;
};

} // end chemkit namespace

#include "graphicsvector-inline.h"

#endif // CHEMKIT_GRAPHICSVECTOR_H
