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

#ifndef CHEMKIT_GRAPHICSMATERIAL_H
#define CHEMKIT_GRAPHICSMATERIAL_H

#include "graphics.h"

namespace chemkit {

class GraphicsMaterialPrivate;

class CHEMKIT_GRAPHICS_EXPORT GraphicsMaterial
{
    public:
        // construction and destruction
        GraphicsMaterial();
        ~GraphicsMaterial();

        // properties
        void setAmbientColor(const QColor &ambient);
        QColor ambientColor() const;
        void setDiffuseColor(const QColor &diffuse);
        QColor diffuseColor() const;
        void setSpecularColor(const QColor &specular);
        QColor specularColor() const;
        void setShininess(int shininess);
        int shininess() const;

    private:
        GraphicsMaterialPrivate* const d;
};

} // end chemkit namespace

#endif // CHEMKIT_GRAPHICSMATERIAL_H
