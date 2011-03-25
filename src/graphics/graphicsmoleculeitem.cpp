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

#include "graphicsmoleculeitem.h"

#include "graphicsscene.h"

#include <chemkit/molecule.h>

namespace chemkit {

// === GraphicsMoleculeItemPrivate ========================================= //
class GraphicsMoleculeItemPrivate
{
    public:
        const Molecule *molecule;
        MoleculeWatcher *watcher;
        GraphicsMoleculeItem::DisplayType displayType;
        float atomRadius;
        float bondRadius;
        float hydrogenScale;
        bool hydrogensVisible;
        bool bondOrderVisible;
        bool atomColoredBonds;
        GraphicsAtomColorMap *colorMap;
        QList<GraphicsAtomItem *> atomItems;
        QList<GraphicsBondItem *> bondItems;
        QList<const Atom *> hiddenAtoms;
};

// === GraphicsMoleculeItem ================================================ //
/// \class GraphicsMoleculeItem graphicsmoleculeitem.h chemkit/graphicsmoleculeitem.h
/// \ingroup chemkit-graphics
/// \brief The GraphicsMoleculeItem visually represents a molecule.
///
/// The GraphicsMoleculeItem class can display molecules in three
/// different ways. The setDisplayType() method is used to switch
/// between them.
///     - \b Ball \b and \b Stick: Each atom is represented by a
///       sphere and each bond by a cylinder.
///     - \b Stick: Every bond is represented by a cylinder and each
///       atom is represented as the cap of a cylider.
///     - \b Space \b Filling: Each atom is represented by a sphere
///       with a radius corresponding to its van der waals radius.
///       This mode is also known as the CPK model.
///
/// \image html molecule-item-types.png

/// \enum GraphicsMoleculeItem::DisplayType
/// Provides names for the different display types:
///     - \c BallAndStick
///     - \c Stick
///     - \c SpaceFilling

// --- Construction and Destruction ---------------------------------------- //
/// Creates a new molecule item to display \p molecule.
GraphicsMoleculeItem::GraphicsMoleculeItem(const Molecule *molecule)
    : QObject(),
      GraphicsItem(MoleculeItem),
      d(new GraphicsMoleculeItemPrivate)
{
    d->molecule = 0;
    d->watcher = new MoleculeWatcher;
    d->atomRadius = 0.5f;
    d->bondRadius = 0.15f;
    d->hydrogenScale = 0.7f;
    d->hydrogensVisible = true;
    d->bondOrderVisible = true;
    d->atomColoredBonds = true;
    d->displayType = BallAndStick;
    d->colorMap = new GraphicsAtomColorMap(GraphicsAtomColorMap::DefaultColorScheme);

    connect(d->watcher, SIGNAL(atomAdded(const chemkit::Atom*)), SLOT(atomAdded(const chemkit::Atom*)));
    connect(d->watcher, SIGNAL(atomRemoved(const chemkit::Atom*)), SLOT(atomRemoved(const chemkit::Atom*)));
    connect(d->watcher, SIGNAL(atomAtomicNumberChanged(const chemkit::Atom*)), SLOT(atomAtomicNumberChanged(const chemkit::Atom*)));
    connect(d->watcher, SIGNAL(atomPositionChanged(const chemkit::Atom*)), SLOT(atomPositionChanged(const chemkit::Atom*)));
    connect(d->watcher, SIGNAL(bondAdded(const chemkit::Bond*)), SLOT(bondAdded(const chemkit::Bond*)));
    connect(d->watcher, SIGNAL(bondRemoved(const chemkit::Bond*)), SLOT(bondRemoved(const chemkit::Bond*)));
    connect(d->watcher, SIGNAL(bondOrderChanged(const chemkit::Bond*)), SLOT(bondOrderChanged(const chemkit::Bond*)));

    setMolecule(molecule);
}

/// Destroys the molecule item.
GraphicsMoleculeItem::~GraphicsMoleculeItem()
{
    qDeleteAll(d->atomItems);
    qDeleteAll(d->bondItems);

    delete d->watcher;
    delete d->colorMap;
}

// --- Properties ---------------------------------------------------------- //
/// Sets the molecule for the item to display to \p molecule.
void GraphicsMoleculeItem::setMolecule(const Molecule *molecule)
{
    d->molecule = molecule;

    qDeleteAll(d->atomItems);
    d->atomItems.clear();
    qDeleteAll(d->bondItems);
    d->bondItems.clear();

    if(molecule){
        foreach(const Atom *atom, molecule->atoms()){
            atomAdded(atom);
        }
        foreach(const Bond *bond, molecule->bonds()){
            bondAdded(bond);
        }
    }

    d->watcher->setMolecule(molecule);
}

/// Returns the molecule that the item displays.
const Molecule* GraphicsMoleculeItem::molecule() const
{
    return d->molecule;
}

/// Set the display type for the molecule.
void GraphicsMoleculeItem::setDisplayType(DisplayType type)
{
    d->displayType = type;

    if(type == BallAndStick){
        setAtomRadius(0.5f);
        setBondRadius(0.15f);
    }
    else if(type == Stick){
        setAtomRadius(0.15f);
        setBondRadius(0.15f);
    }
    else if(type == SpaceFilling){
        foreach(GraphicsAtomItem *item, d->atomItems){
            item->setRadius(item->atom()->vanDerWaalsRadius());
        }
    }

    update();
}

/// Returns the current display type used for the molecule.
GraphicsMoleculeItem::DisplayType GraphicsMoleculeItem::displayType() const
{
    return d->displayType;
}

/// Sets the radius of the spheres used for displaying the atoms to
/// \p radius.
void GraphicsMoleculeItem::setAtomRadius(float radius)
{
    d->atomRadius = radius;

    foreach(GraphicsAtomItem *item, d->atomItems){
        if(d->displayType == BallAndStick && item->atom()->isTerminalHydrogen())
            item->setRadius(radius * d->hydrogenScale);
        else
            item->setRadius(radius);
    }

    foreach(GraphicsBondItem *item, d->bondItems){
        if(d->displayType == BallAndStick && item->bond()->isTerminal() && item->bond()->contains(Atom::Hydrogen))
            item->setMaximumRadius(radius * d->hydrogenScale);
        else
            item->setMaximumRadius(radius);
    }
}

/// Returns the radius of the spheres used for displaying atoms.
float GraphicsMoleculeItem::atomRadius() const
{
    return d->atomRadius;
}

/// Sets the radius of the cylinders used for displaing bonds to
/// \p radius.
void GraphicsMoleculeItem::setBondRadius(float radius)
{
    d->bondRadius = radius;

    foreach(GraphicsBondItem *item, d->bondItems){
        item->setRadius(radius);
    }
}

/// Returns the radius of the cylinders used for displaying bonds.
float GraphicsMoleculeItem::bondRadius() const
{
    return d->bondRadius;
}

/// Sets whether or not to show terminal hydrogen atoms.
void GraphicsMoleculeItem::setHydrogensVisible(bool visible)
{
    d->hydrogensVisible = visible;

    foreach(GraphicsAtomItem *item, d->atomItems){
        if(item->atom()->isTerminalHydrogen()){
            item->setVisible(visible);
        }
    }
    foreach(GraphicsBondItem *item, d->bondItems){
        if(item->bond()->contains(Atom::Hydrogen)){
            item->setVisible(visible);
        }
    }

    update();
}

/// Returns \c true if terminal hydrogen atoms are being shown.
bool GraphicsMoleculeItem::hydrogensVisible() const
{
    return d->hydrogensVisible;
}

void GraphicsMoleculeItem::setHydrogenScale(float scale)
{
    d->hydrogenScale = scale;
    update();
}

float GraphicsMoleculeItem::hydrogenScale() const
{
    return d->hydrogenScale;
}

void GraphicsMoleculeItem::setBondOrderVisible(bool showBondOrder)
{
    foreach(GraphicsBondItem *item, d->bondItems){
        item->setBondOrderVisible(showBondOrder);
    }

    d->bondOrderVisible = showBondOrder;
    update();
}

bool GraphicsMoleculeItem::bondOrderVisible() const
{
    return d->bondOrderVisible;
}

void GraphicsMoleculeItem::setAtomColoredBonds(bool atomColoredBonds)
{
    foreach(GraphicsBondItem *item, d->bondItems){
        item->setAtomColored(atomColoredBonds);
    }

    d->atomColoredBonds = atomColoredBonds;
    update();
}

bool GraphicsMoleculeItem::atomColoredBonds() const
{
    return d->atomColoredBonds;
}

void GraphicsMoleculeItem::setAtomVisible(const Atom *atom, bool visible)
{
    if(visible){
        d->hiddenAtoms.removeAll(atom);
    }
    else{
        d->hiddenAtoms.append(atom);
    }
}

bool GraphicsMoleculeItem::atomVisible(const Atom *atom) const
{
    return !d->hiddenAtoms.contains(atom);
}

/// Sets the color map for the molecule item to \p colorMap.
void GraphicsMoleculeItem::setAtomColorMap(GraphicsAtomColorMap *colorMap)
{
    delete d->colorMap;

    d->colorMap = colorMap;
}

/// Returns the color map for the molecule item.
GraphicsAtomColorMap* GraphicsMoleculeItem::atomColorMap() const
{
    return d->colorMap;
}

// --- Items --------------------------------------------------------------- //
GraphicsAtomItem* GraphicsMoleculeItem::atomItem(const Atom *atom) const
{
    foreach(GraphicsAtomItem *item, d->atomItems){
        if(item->atom() == atom){
            return item;
        }
    }

    return 0;
}

GraphicsBondItem* GraphicsMoleculeItem::bondItem(const Bond *bond) const
{
    foreach(GraphicsBondItem *item, d->bondItems){
        if(item->bond() == bond){
            return item;
        }
    }

    return 0;
}

// --- Painting ------------------------------------------------------------ //
void GraphicsMoleculeItem::paint(GraphicsPainter *painter)
{
    Q_UNUSED(painter);
}

// --- Events -------------------------------------------------------------- //
void GraphicsMoleculeItem::itemChanged(ItemChange change)
{
    if(change == ItemSceneChanged){
        foreach(GraphicsItem *item, d->atomItems){
            if(item->scene())
                item->scene()->removeItem(item);
            if(scene())
                scene()->addItem(item);
        }
        foreach(GraphicsItem *item, d->bondItems){
            if(item->scene())
                item->scene()->removeItem(item);
            if(scene())
                scene()->addItem(item);
        }
    }
    else if(change == ItemVisiblityChanged){
        foreach(GraphicsItem *item, d->atomItems){
            item->setVisible(isVisible());
        }
        foreach(GraphicsItem *item, d->bondItems){
            item->setVisible(isVisible());
        }
    }
}

// --- Slots --------------------------------------------------------------- //
void GraphicsMoleculeItem::atomAdded(const chemkit::Atom *atom)
{
    float radius;
    if(d->displayType == SpaceFilling){
        radius = atom->vanDerWaalsRadius();
    }
    else{
        radius = d->atomRadius;

        if(d->displayType == BallAndStick && atom->isTerminalHydrogen()){
            radius *= d->hydrogenScale;
        }
    }

    GraphicsAtomItem *item = new GraphicsAtomItem(const_cast<Atom *>(atom), radius);
    item->setVisible(isVisible());
    item->setColor(d->colorMap->color(atom));

    if(scene()){
        scene()->addItem(item);
    }

    d->atomItems.append(item);
}

void GraphicsMoleculeItem::atomRemoved(const chemkit::Atom *atom)
{
    foreach(GraphicsAtomItem *item, d->atomItems){
        if(item->atom() == atom){
            d->atomItems.removeOne(item);

            if(scene()){
                scene()->removeItem(item);
            }

            delete item;

            break;
        }
    }
}

void GraphicsMoleculeItem::atomAtomicNumberChanged(const chemkit::Atom *atom)
{
    GraphicsAtomItem *item = atomItem(atom);
    item->setColor(d->colorMap->color(atom));

    if(d->displayType == SpaceFilling){
        item->setRadius(atom->vanDerWaalsRadius());
    }
    else if(d->displayType == BallAndStick && atom->isTerminalHydrogen()){
        item->setRadius(d->atomRadius * d->hydrogenScale);
    }
    else{
        item->setRadius(d->atomRadius);
    }

    foreach(const Bond *bond, atom->bonds()){
        GraphicsBondItem *item = bondItem(bond);
        item->setAtomColors(d->colorMap->color(bond->atom1()),
                            d->colorMap->color(bond->atom2()));
    }

    update();
}

void GraphicsMoleculeItem::atomPositionChanged(const chemkit::Atom *atom)
{
    foreach(GraphicsAtomItem *item, d->atomItems){
        if(item->atom() == atom){
            item->setAtom(atom);
            item->update();
        }
    }

    foreach(GraphicsBondItem *item, d->bondItems){
        if(item->bond()->contains(atom)){
            item->update();
        }
    }
}

void GraphicsMoleculeItem::bondAdded(const chemkit::Bond *bond)
{
    GraphicsBondItem *item = new GraphicsBondItem(const_cast<Bond *>(bond));
    item->setVisible(isVisible());

    item->setAtomColors(d->colorMap->color(bond->atom1()),
                        d->colorMap->color(bond->atom2()));

    if(scene()){
        scene()->addItem(item);
    }

    d->bondItems.append(item);

    if(d->displayType == BallAndStick && bond->isTerminal() && bond->contains(Atom::Hydrogen)){
        const Atom *hydrogen = bond->atom1()->isTerminalHydrogen() ? bond->atom1() : bond->atom2();
        GraphicsAtomItem *item = atomItem(hydrogen);
        if(item){
            item->setRadius(d->atomRadius * d->hydrogenScale);
        }
    }

    update();
}

void GraphicsMoleculeItem::bondRemoved(const chemkit::Bond *bond)
{
    foreach(GraphicsBondItem *item, d->bondItems){
        if(item->bond() == bond){
            d->bondItems.removeOne(item);

            if(scene())
                scene()->removeItem(item);
            delete item;
            break;
        }
    }

    update();
}

void GraphicsMoleculeItem::bondOrderChanged(const chemkit::Bond *bond)
{
    Q_UNUSED(bond);

    update();
}

} // end chemkit namespace
