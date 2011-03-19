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

#ifndef CHEMKIT_VECTOR3G_H
#define CHEMKIT_VECTOR3G_H

#include "graphics.h"

#include <chemkit/vector3.h>
#include <chemkit/genericvector.h>

namespace chemkit {

class CHEMKIT_GRAPHICS_EXPORT Vector3g : public GenericVector<GraphicsFloat>
{
    public:
        // construction and destruction
        Vector3g();
        Vector3g(GraphicsFloat x, GraphicsFloat y, GraphicsFloat z);
        Vector3g(const GenericVector<float> &vector);
        Vector3g(const GenericVector<double> &vector);
        Vector3g(const StaticVector<float, 3> &vector);
        Vector3g(const StaticVector<double, 3> &vector);
};

} // end chemkit namespace

#include "vector3g-inline.h"

#endif // CHEMKIT_VECTOR3G_H
