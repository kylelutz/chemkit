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

#include "graphicsmaterial.h"

namespace chemkit {

// === GraphicsMaterialPrivate ============================================= //
class GraphicsMaterialPrivate
{
    public:
        int shininess;
        QColor ambientColor;
        QColor diffuseColor;
        QColor specularColor;
};

// === GraphicsMaterial ==================================================== //
/// \class GraphicsMaterial graphicsmaterial.h chemkit/graphicsmaterial.h
/// \ingroup chemkit-graphics
/// \internal
/// \brief The GraphicsMaterial class represents a material.

// --- Construction and Destruction ---------------------------------------- //
/// Creates a new graphics material.
GraphicsMaterial::GraphicsMaterial()
    : d(new GraphicsMaterialPrivate)
{
    d->shininess = 15;
    d->specularColor = QColor::fromRgbF(0.3f, 0.3f, 0.3f, 1.0f);
}

/// Destroys the graphics material object.
GraphicsMaterial::~GraphicsMaterial()
{
    delete d;
}

// --- Properties ---------------------------------------------------------- //
/// Sets the ambient color to \p color.
void GraphicsMaterial::setAmbientColor(const QColor &color)
{
    d->ambientColor = color;
}

/// Returns the ambient color.
QColor GraphicsMaterial::ambientColor() const
{
    return d->ambientColor;
}

/// Sets the diffuse color to \p color.
void GraphicsMaterial::setDiffuseColor(const QColor &color)
{
    d->diffuseColor = color;
}

/// Returns the diffuse color.
QColor GraphicsMaterial::diffuseColor() const
{
    return d->diffuseColor;
}

/// Sets the specular color to \p color.
void GraphicsMaterial::setSpecularColor(const QColor &color)
{
    d->specularColor = color;
}

/// Returns the specular color.
QColor GraphicsMaterial::specularColor() const
{
    return d->specularColor;
}

/// Sets the shininess to \p shininess.
void GraphicsMaterial::setShininess(int shininess)
{
    d->shininess = shininess;
}

/// Returns the shininess.
int GraphicsMaterial::shininess() const
{
    return d->shininess;
}

} // end chemkit namespace
