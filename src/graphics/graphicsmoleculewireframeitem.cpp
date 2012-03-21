/******************************************************************************
**
** Copyright (C) 2009-2012 Kyle Lutz <kyle.r.lutz@gmail.com>
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

#include "graphicsmoleculewireframeitem.h"

#include <boost/make_shared.hpp>

#include <chemkit/atom.h>
#include <chemkit/bond.h>
#include <chemkit/foreach.h>
#include <chemkit/molecule.h>
#include <chemkit/atomcolormap.h>

#include "graphicspainter.h"
#include "graphicsvertexbuffer.h"

namespace chemkit {

// === GraphicsMoleculeWireframeItemPrivate ================================ //
class GraphicsMoleculeWireframeItemPrivate
{
public:
    const Molecule *molecule;
    boost::shared_ptr<AtomColorMap> colorMap;
    bool hydrogensVisibile;
};

// === GraphicsMoleculeWireframeItem ======================================= //
/// \class GraphicsMoleculeWireframeItem graphicsmoleculewireframeitem.h chemkit/graphicsmoleculewireframeitem.h
/// \ingroup chemkit-graphics
/// \brief The GraphicsMoleculeWireframeItem shows a molecule as a
///        wireframe.

// --- Construction and Destruction ---------------------------------------- //
/// Creates a new molecule wireframe item for \p molecule.
GraphicsMoleculeWireframeItem::GraphicsMoleculeWireframeItem(const Molecule *molecule)
    : d(new GraphicsMoleculeWireframeItemPrivate)
{
    d->molecule = molecule;
    d->colorMap = boost::make_shared<AtomColorMap>(AtomColorMap::DefaultColorScheme);
    d->hydrogensVisibile = true;
}

/// Destroys the molecule wireframe item.
GraphicsMoleculeWireframeItem::~GraphicsMoleculeWireframeItem()
{
    delete d;
}

// --- Properties ---------------------------------------------------------- //
/// Sets the molecule to display to \p molecule.
void GraphicsMoleculeWireframeItem::setMolecule(const Molecule *molecule)
{
    d->molecule = molecule;
}

/// Returns the molecule currently being displayed.
const Molecule* GraphicsMoleculeWireframeItem::molecule() const
{
    return d->molecule;
}

/// Sets the color map to \p colorMap.
void GraphicsMoleculeWireframeItem::setColorMap(const boost::shared_ptr<AtomColorMap> &colorMap)
{
    d->colorMap = colorMap;
}

/// Returns the current color map.
boost::shared_ptr<AtomColorMap> GraphicsMoleculeWireframeItem::colorMap() const
{
    return d->colorMap;
}

/// Enables/disables displaying terminal hydrogen atoms.
void GraphicsMoleculeWireframeItem::setHydrogensVisibible(bool visibile)
{
    d->hydrogensVisibile = visibile;
}

/// Returns \c true if terminal hydrogens are being displayed.
bool GraphicsMoleculeWireframeItem::hydrogensVisibile()
{
    return d->hydrogensVisibile;
}

// --- Drawing ------------------------------------------------------------- //
void GraphicsMoleculeWireframeItem::paint(GraphicsPainter *painter)
{
    if(!d->molecule){
        return;
    }

    QVector<Point3f> vertices;
    QVector<unsigned short> indices;
    QVector<QColor> colors;

    foreach(const Atom *atom, d->molecule->atoms()){
        vertices.append(atom->position().cast<float>());
        colors.append(d->colorMap->color(atom));
    }

    foreach(const Bond *bond, d->molecule->bonds()){
        if(!d->hydrogensVisibile &&
           bond->isTerminal() &&
           bond->contains(Atom::Hydrogen)){
            continue;
        }

        const Atom *a = bond->atom1();
        const Atom *b = bond->atom2();

        QColor colorA = colors[a->index()];
        QColor colorB = colors[b->index()];

        if(colorA == colorB){
            indices.append(a->index());
            indices.append(b->index());
        }
        else{
            Point3f center = bond->center().cast<float>();
            size_t centerIndex = vertices.size();
            vertices.append(center);
            colors.append(colorA);

            // a -> center
            indices.append(a->index());
            indices.append(centerIndex);

            // center -> b
            indices.append(centerIndex);
            indices.append(b->index());
        }
    }

    GraphicsVertexBuffer buffer;

    buffer.setVertices(vertices);
    buffer.setIndices(indices);
    buffer.setColors(colors);

    glDisable(GL_LIGHTING);
    glShadeModel(GL_FLAT);

    QGLShaderProgram program;
    program.addShaderFromSourceFile(QGLShader::Vertex, ":/shaders/flat.vert");
    program.addShaderFromSourceFile(QGLShader::Fragment, ":/shaders/flat.frag");
    program.link();
    program.bind();

    painter->draw(&buffer, GraphicsPainter::Lines);

    program.release();

    glShadeModel(GL_SMOOTH);
    glEnable(GL_LIGHTING);
}

} // end chemkit namespace
