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

#include "graphicsmolecularsurfaceitem.h"

#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <boost/make_shared.hpp>

#include <chemkit/atom.h>
#include <chemkit/foreach.h>
#include <chemkit/geometry.h>
#include <chemkit/molecule.h>
#include <chemkit/alphashape.h>

#include "graphicsscene.h"
#include "graphicssphere.h"
#include "graphicspainter.h"
#include "graphicsmaterial.h"
#include "graphicsmoleculeitem.h"
#include "graphicsvertexbuffer.h"

namespace chemkit {

namespace {

// Returns the electrostatic potential at the given position calculated
// from the partial charges and positions of the atoms in the molecule.
float electrostaticPotential(const Molecule *molecule, const Point3f &position)
{
    float esp = 0;

    const float pi = chemkit::constants::Pi;
    const float e0 = 1.0f;

    foreach(const Atom *atom, molecule->atoms()){
        float q = atom->partialCharge();
        float r = (position - atom->position().cast<float>()).norm();

        esp += (1.0f / (4.0f * pi * e0)) * (q / r);
    }

    return esp;
}

// Returns a color interpolated at the given value between the colors a
// (starting at value av) and b (starting at value bv).
QColor interpolate(const QColor &a, const QColor &b, float av, float bv, float value)
{
    float v = (value - av) / (bv - av);

    return QColor::fromRgbF(a.redF() + (b.redF() - a.redF()) * v,
                            a.greenF() + (b.greenF() - a.greenF()) * v,
                            a.blueF() + (b.blueF() - a.blueF()) * v);
}

// Returns the color associated with the electrostatic potential.
QColor electrostaticPotentialColor(float esp)
{
    QColor red = QColor::fromRgb(255, 0, 0);
    QColor orange = QColor::fromRgb(255, 127, 0);
    QColor yellow = QColor::fromRgb(255, 255, 0);
    QColor green = QColor::fromRgb(0, 255, 0);
    QColor blue = QColor::fromRgb(0, 0, 255);

    // color ranges are hard coded (for now)
    float redStart = -0.0075f;
    float orangeStart = -0.0035f;
    float yellowStart = 0.0f;
    float greenStart = 0.0015f;
    float blueStart = 0.0045f;

    if(esp < redStart){
        return red;
    }
    else if(esp < orangeStart){
        return interpolate(red, orange, redStart, orangeStart, esp);
    }
    else if(esp < yellowStart){
        return interpolate(orange, yellow, orangeStart, yellowStart, esp);
    }
    else if(esp < greenStart){
        return interpolate(yellow, green, yellowStart, greenStart, esp);
    }
    else if(esp < blueStart){
        return interpolate(green, blue, greenStart, blueStart, esp);
    }
    else{
        return blue;
    }
}

} // end anonymous namespace

// === GraphicsMolecularSurfaceItemPrivate ================================= //
class GraphicsMolecularSurfaceItemPrivate
{
public:
    MolecularSurface *surface;
    QColor color;
    boost::shared_ptr<AtomColorMap> colorMap;
    GraphicsMolecularSurfaceItem::ColorMode colorMode;
};

// === GraphicsMolecularSurfaceItem ======================================== //
/// \class GraphicsMolecularSurfaceItem graphicsmolecularsurfaceitem.h chemkit/graphicsmolecularsurfaceitem.h
/// \ingroup chemkit-graphics
/// \brief The GraphicsMolecularSurfaceItem displays a molecular
///        surface.
///
/// The image below shows a molecular surface item on top of a
/// guanine molecule:
/// \image html molecular-surface-item.png
///
/// \see MolecularSurface

/// \enum GraphicsMolecularSurfaceItem::ColorMode
/// Provides different modes for coloring the surface.
///     - \c SolidColor
///     - \c AtomColor
///     - \c ElectrostaticPotential

// --- Construction and Destruction ---------------------------------------- //
/// Creates a new molecular surface item for \p molecule.
GraphicsMolecularSurfaceItem::GraphicsMolecularSurfaceItem(const Molecule *molecule)
    : GraphicsItem(),
      d(new GraphicsMolecularSurfaceItemPrivate)
{
    d->surface = new MolecularSurface(molecule, MolecularSurface::SolventExcluded);
    d->color = Qt::red;
    d->colorMode = AtomColor;
    d->colorMap = boost::make_shared<AtomColorMap>(AtomColorMap::DefaultColorScheme);

    setMolecule(molecule);
}

/// Creates a new molecular surface item for \p surface.
GraphicsMolecularSurfaceItem::GraphicsMolecularSurfaceItem(const MolecularSurface *surface)
    : GraphicsItem(),
      d(new GraphicsMolecularSurfaceItemPrivate)
{
    d->surface = new MolecularSurface(surface->molecule(), MolecularSurface::SolventExcluded);
    d->color = Qt::red;
    d->colorMode = AtomColor;
    d->colorMap = boost::make_shared<AtomColorMap>(AtomColorMap::DefaultColorScheme);

    setSurface(surface);
}

/// Destroys the molecular surface item.
GraphicsMolecularSurfaceItem::~GraphicsMolecularSurfaceItem()
{
    delete d->surface;
    delete d;
}

// --- Properties ---------------------------------------------------------- //
/// Sets the surface to display to \p surface.
void GraphicsMolecularSurfaceItem::setSurface(const MolecularSurface *surface)
{
    if(surface){
        d->surface->setMolecule(surface->molecule());
        d->surface->setSurfaceType(surface->surfaceType());
        d->surface->setProbeRadius(surface->probeRadius());
    }
    else{
        d->surface->setMolecule(0);
    }
}

/// Returns the surface that the item is displaying.
const MolecularSurface* GraphicsMolecularSurfaceItem::surface() const
{
    return d->surface;
}

/// Sets the molecule for the surface to \p molecule.
void GraphicsMolecularSurfaceItem::setMolecule(const Molecule *molecule)
{
    d->surface->setMolecule(molecule);
}

/// Returns the molecule for the surface.
const Molecule* GraphicsMolecularSurfaceItem::molecule() const
{
    return d->surface->molecule();
}

/// Sets the surface type to \p type.
void GraphicsMolecularSurfaceItem::setSurfaceType(MolecularSurface::SurfaceType type)
{
    d->surface->setSurfaceType(type);
}

/// Returns the surface type.
MolecularSurface::SurfaceType GraphicsMolecularSurfaceItem::surfaceType() const
{
    return d->surface->surfaceType();
}

/// Sets the probe radius for the surface to \p radius.
void GraphicsMolecularSurfaceItem::setProbeRadius(float radius)
{
    d->surface->setProbeRadius(radius);
}

/// Returns the probe radius for the surface.
float GraphicsMolecularSurfaceItem::probeRadius() const
{
    return d->surface->probeRadius();
}

/// Sets the color for the surface to \p color.
void GraphicsMolecularSurfaceItem::setColor(const QColor &color)
{
    d->color = color;
}

/// Returns the color for the surface.
QColor GraphicsMolecularSurfaceItem::color() const
{
    return d->color;
}

/// Sets the color mode for the surface item to \p mode.
void GraphicsMolecularSurfaceItem::setColorMode(ColorMode mode)
{
    d->colorMode = mode;
}

/// Returns the color mode for the surface item.
GraphicsMolecularSurfaceItem::ColorMode GraphicsMolecularSurfaceItem::colorMode() const
{
    return d->colorMode;
}

/// Sets the color map for the surface item to \p colorMap.
void GraphicsMolecularSurfaceItem::setColorMap(const boost::shared_ptr<AtomColorMap> &colorMap)
{
    d->colorMap = colorMap;
}

/// Returns the color map for the surface item.
boost::shared_ptr<AtomColorMap> GraphicsMolecularSurfaceItem::colorMap() const
{
    return d->colorMap;
}

// --- Drawing ------------------------------------------------------------- //
void GraphicsMolecularSurfaceItem::paint(GraphicsPainter *painter)
{
    if(!d->surface || !d->surface->molecule()){
        // nothing to render
        return;
    }

    // create map of sphere radii -> sphere vertex buffers
    std::map<Real, boost::shared_ptr<GraphicsVertexBuffer> > sphereCache;

    foreach(const Atom *atom, d->surface->molecule()->atoms()){
        if(sphereCache.find(atom->vanDerWaalsRadius()) == sphereCache.end()){
            float radius = atom->vanDerWaalsRadius();
            if(d->surface->surfaceType() != MolecularSurface::VanDerWaals){
                radius += d->surface->probeRadius();
            }

            GraphicsSphere sphere(radius);

            sphereCache[atom->vanDerWaalsRadius()] =
                boost::shared_ptr<GraphicsVertexBuffer>(sphere.tesselate());
        }
    }

    // draw opaque
    if(isOpaque()){
        foreach(const Atom *atom, d->surface->molecule()->atoms()){
            const boost::shared_ptr<GraphicsVertexBuffer> &buffer =
                sphereCache.find(atom->vanDerWaalsRadius())->second;

            if(d->colorMode == AtomColor){
                painter->setColor(d->colorMap->color(atom));
            }
            else if(d->colorMode == SolidColor){
                painter->setColor(color());
            }
            else if(d->colorMode == ElectrostaticPotential){
                // function that returns electrostatic potential at a point in space
                boost::function<float (const Point3f&)> espFunction =
                    boost::bind(electrostaticPotential, d->surface->molecule(), _1);

                // calculate color for each vertex
                QVector<QColor> colors;
                foreach(const Point3f &vertex, buffer->vertices()){
                    float esp = espFunction(atom->position().cast<float>() + vertex);
                    QColor color = electrostaticPotentialColor(esp);
                    colors.append(color);
                }

                buffer->setColors(colors);
            }

            glPushMatrix();
            glTranslated(atom->x(), atom->y(), atom->z());
            painter->draw(buffer.get());
            glPopMatrix();
        }
    }
    // draw translucent
    else{
        // disable writes to color buffer
        glColorMask(false, false, false, false);

        // render each sphere to fill depth buffer
        foreach(const Atom *atom, d->surface->molecule()->atoms()){
            const boost::shared_ptr<GraphicsVertexBuffer> &buffer =
                sphereCache.find(atom->vanDerWaalsRadius())->second;

            glPushMatrix();
            glTranslated(atom->x(), atom->y(), atom->z());
            painter->draw(buffer.get());
            glPopMatrix();
        }

        // re-enabled writes to color buffer
        glColorMask(true, true, true, true);

        // only select the closest fragments
        glDepthFunc(GL_EQUAL);

        // disable writes to the depth buffer
        glDepthMask(false);

        // actually render each sphere
        foreach(const Atom *atom, d->surface->molecule()->atoms()){
            const boost::shared_ptr<GraphicsVertexBuffer> &buffer =
                sphereCache.find(atom->vanDerWaalsRadius())->second;

            if(d->colorMode == AtomColor){
                QColor color = d->colorMap->color(atom);
                color.setAlphaF(opacity());
                painter->setColor(color);
            }
            else if(d->colorMode == SolidColor){
                QColor color = this->color();
                color.setAlphaF(opacity());
                painter->setColor(color);
            }
            else if(d->colorMode == ElectrostaticPotential){
                // function that returns electrostatic potential at a point in space
                boost::function<float (const Point3f&)> espFunction =
                    boost::bind(electrostaticPotential, d->surface->molecule(), _1);

                // calculate color for each vertex
                QVector<QColor> colors;
                foreach(const Point3f &vertex, buffer->vertices()){
                    float esp = espFunction(atom->position().cast<float>() + vertex);
                    QColor color = electrostaticPotentialColor(esp);
                    color.setAlphaF(opacity());
                    colors.append(color);
                }

                buffer->setColors(colors);
            }

            glPushMatrix();
            glTranslated(atom->x(), atom->y(), atom->z());
            painter->draw(buffer.get());
            glPopMatrix();
        }

        // reset depth buffer
        glDepthFunc(GL_LESS);
        glDepthMask(true);
    }
}

} // end chemkit namespace
