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

#include "graphicsbonditem.h"

#include <chemkit/atom.h>
#include <chemkit/bond.h>
#include <chemkit/point3.h>
#include <chemkit/vector3.h>
#include <chemkit/geometry.h>

#include "graphicsray.h"
#include "graphicspainter.h"
#include "graphicsmoleculeitem.h"

namespace chemkit {

// === GraphicsBondItemPrivate ============================================= //
class GraphicsBondItemPrivate
{
public:
    const Bond *bond;
    float radius;
    float maximumRadius;
    Vector3f normal;
    bool bondOrderVisible;
    bool atomColored;
    QColor color;
    QPair<QColor, QColor> atomColors;
};

// === GraphicsBondItem ==================================================== //
/// \class GraphicsBondItem graphicsbonditem.h chemkit/graphicsbonditem.h
/// \ingroup chemkit-graphics
/// \brief The GraphicsBondItem visually represents a chemical bond.
///
/// The GraphicsBondItem class displays a single bond using a
/// cylinder.

// --- Construction and Destruction ---------------------------------------- //
/// Creates a new bond item for \p bond with a cylinder of the given
/// \p radius.
GraphicsBondItem::GraphicsBondItem(const Bond *bond, float radius)
    : GraphicsItem(BondItem),
      d(new GraphicsBondItemPrivate)
{
    d->bond = bond;
    d->radius = radius;
    d->maximumRadius = 0.5;
    d->normal = -Vector3f::UnitZ();
    d->bondOrderVisible = true;
    d->atomColored = true;
    d->color = Qt::darkGray;

    setBond(bond);
}

/// Destoys the bond item.
GraphicsBondItem::~GraphicsBondItem()
{
    delete d;
}

// --- Properties ---------------------------------------------------------- //
/// Sets the bond item that is displayed to \p bond.
void GraphicsBondItem::setBond(const Bond *bond)
{
    d->bond = bond;
}

/// Returns the bond being displayed by the item or \c 0 if no bond
/// is being displayed.
const Bond* GraphicsBondItem::bond() const
{
    return d->bond;
}

/// Sets the radius of the cylinder to \p radius.
void GraphicsBondItem::setRadius(float radius)
{
    d->radius = radius;
}

/// Returns the radius of the cylinder.
float GraphicsBondItem::radius() const
{
    return d->radius;
}

void GraphicsBondItem::setMaximumRadius(float radius)
{
    d->maximumRadius = radius;
}

float GraphicsBondItem::maximumRadius() const
{
    return d->maximumRadius;
}

void GraphicsBondItem::setNormal(const Vector3f &normal)
{
    d->normal = normal;
}

Vector3f GraphicsBondItem::normal() const
{
    return d->normal;
}

void GraphicsBondItem::setAtomColored(bool atomColored)
{
    d->atomColored = atomColored;
}

bool GraphicsBondItem::atomColored() const
{
    return d->atomColored;
}

void GraphicsBondItem::setBondOrderVisible(bool visible)
{
    d->bondOrderVisible = visible;
}

bool GraphicsBondItem::bondOrderVisible() const
{
    return d->bondOrderVisible;
}

/// Sets the color for the bond to \p color.
void GraphicsBondItem::setColor(const QColor &color)
{
    d->color = color;
}

/// Returns the color for the bond.
QColor GraphicsBondItem::color() const
{
    return d->color;
}

/// Sets the color for the bond to \p a and \p b.
void GraphicsBondItem::setAtomColors(const QColor &a, const QColor &b)
{
    d->atomColors = qMakePair(a, b);
}

/// Returns the atom colors for the bond.
QPair<QColor, QColor> GraphicsBondItem::atomColors()
{
    return d->atomColors;
}

// --- Intersection -------------------------------------------------------- //
bool GraphicsBondItem::intersects(const GraphicsRay &ray, float *distance) const
{
    float intersectionRadius = qMin(d->radius, d->maximumRadius);

    return ray.intersectsCylinder(d->bond->atom1()->position().cast<float>(),
                                  d->bond->atom2()->position().cast<float>(),
                                  intersectionRadius,
                                  distance);
}

// --- Drawing ------------------------------------------------------------- //
void GraphicsBondItem::paint(GraphicsPainter *painter)
{
    const Atom *atom1 = d->bond->atom1();
    const Atom *atom2 = d->bond->atom2();
    int bondOrder = d->bond->order();

    // draw multiple cylinders
    if(d->bondOrderVisible && bondOrder > 1){
        // radius of each cylinder, clamped to be under
        // the total maximum radius if necessary
        float radius = d->radius;
        int radiiCount = 2 * bondOrder - 1;
        if(radius * radiiCount > d->maximumRadius)
            radius = d->maximumRadius / radiiCount;

        // distance between the center's of each adjacent cylinder
        float offset = 3 * radius;

        // distance from the atom's positions to the
        // center of the topmost cylinder
        float initialOffset = 1.5 * (bondOrder - 1) * radius;

        // a vector pointing to the right (with normal pointing up)
        Vector3f right = (atom2->position().cast<float>() - atom1->position().cast<float>()).cross(d->normal);

        // positions for the first cylinder
        Point3f a = atom1->position().cast<float>() + (right.normalized() * -initialOffset);
        Point3f b = atom2->position().cast<float>() + (right.normalized() * -initialOffset);

        // draw each cylinder
        for(int i = 0; i < bondOrder; i++){
            if(d->atomColored){
                if(d->atomColors.first == d->atomColors.second){
                    painter->setColor(d->atomColors.first);
                    painter->drawCylinder(a, b, radius);
                }
                else{
                    Point3f midpoint = chemkit::geometry::midpoint(a.cast<Real>(), b.cast<Real>()).cast<float>();
                    painter->setColor(d->atomColors.first);
                    painter->drawCylinder(a, midpoint, radius);
                    painter->setColor(d->atomColors.second);
                    painter->drawCylinder(midpoint, b, radius);
                }
            }
            else{
                painter->setColor(d->color);
                painter->drawCylinder(a, b, radius);
            }

            // move the positions for the next cylinder
            a += right.normalized() * offset;
            b += right.normalized() * offset;
        }
    }

    // draw a single cylinder
    else{
        float radius = qMin(d->radius, d->maximumRadius);

        if(d->atomColored){
            if(d->atomColors.first == d->atomColors.second){
                painter->setColor(d->atomColors.first);
                painter->drawCylinder(atom1->position().cast<float>(), atom2->position().cast<float>(), radius);
            }
            else{
                Point3f midpoint = chemkit::geometry::midpoint(atom1->position(), atom2->position()).cast<float>();
                painter->setColor(d->atomColors.first);
                painter->drawCylinder(atom1->position().cast<float>(), midpoint, radius);
                painter->setColor(d->atomColors.second);
                painter->drawCylinder(midpoint, atom2->position().cast<float>(), radius);
            }
        }
        else{
            painter->setColor(d->color);
            painter->drawCylinder(atom1->position().cast<float>(), atom2->position().cast<float>(), radius);
        }
    }
}

} // end chemkit namespace
