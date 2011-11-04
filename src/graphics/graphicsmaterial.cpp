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
