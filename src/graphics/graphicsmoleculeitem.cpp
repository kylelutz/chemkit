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
        GraphicsFloat atomRadius;
        GraphicsFloat bondRadius;
        GraphicsFloat hydrogenScale;
        bool hydrogensVisible;
        bool bondOrderVisible;
        bool atomColoredBonds;
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
void GraphicsMoleculeItem::setAtomRadius(GraphicsFloat radius)
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
GraphicsFloat GraphicsMoleculeItem::atomRadius() const
{
    return d->atomRadius;
}

/// Sets the radius of the cylinders used for displaing bonds to
/// \p radius.
void GraphicsMoleculeItem::setBondRadius(GraphicsFloat radius)
{
    d->bondRadius = radius;

    foreach(GraphicsBondItem *item, d->bondItems){
        item->setRadius(radius);
    }
}

/// Returns the radius of the cylinders used for displaying bonds.
GraphicsFloat GraphicsMoleculeItem::bondRadius() const
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

void GraphicsMoleculeItem::setHydrogenScale(GraphicsFloat scale)
{
    d->hydrogenScale = scale;
    update();
}

GraphicsFloat GraphicsMoleculeItem::hydrogenScale() const
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

// --- Items --------------------------------------------------------------- //
GraphicsAtomItem* GraphicsMoleculeItem::atomItem(const Atom *atom)
{
    foreach(GraphicsAtomItem *item, d->atomItems){
        if(item->atom() == atom){
            return item;
        }
    }

    return 0;
}

GraphicsBondItem* GraphicsMoleculeItem::bondItem(const Bond *bond)
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

// --- Static Methods ------------------------------------------------------ //
QColor GraphicsMoleculeItem::atomColor(int atomicNumber)
{
    switch(atomicNumber){
        case 1: return QColor(255, 255, 255);
        case 2: return QColor(217, 255, 255);
        case 3: return QColor(204, 128, 255);
        case 4: return QColor(194, 255, 0);
        case 5: return QColor(255, 181, 181);
        case 6: return QColor(80, 80, 80);
        case 7: return QColor(48, 80, 248);
        case 8: return QColor(255, 13, 13);
        case 9: return QColor(144, 224, 80);
        case 10: return QColor(179, 227, 245);
        case 11: return QColor(171, 92, 242);
        case 12: return QColor(138, 255, 0);
        case 13: return QColor(191, 166, 166);
        case 14: return QColor(240, 200, 160);
        case 15: return QColor(255, 128, 0);
        case 16: return QColor(255, 255, 48);
        case 17: return QColor(31, 240, 31);
        case 18: return QColor(128, 209, 227);
        case 19: return QColor(143, 64, 212);
        case 20: return QColor(61, 255, 0);
        case 35: return QColor(166, 41, 41);
        case 53: return QColor(148, 0, 148);
        default: return QColor(255, 20, 147);
    }
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

QColor GraphicsMoleculeItem::atomColor(const Atom *atom)
{
    return atomColor(atom->atomicNumber());
}

// --- Slots --------------------------------------------------------------- //
void GraphicsMoleculeItem::atomAdded(const chemkit::Atom *atom)
{
    GraphicsFloat radius;
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

    if(d->displayType == SpaceFilling){
        item->setRadius(atom->vanDerWaalsRadius());
    }
    else if(d->displayType == BallAndStick && atom->isTerminalHydrogen()){
        item->setRadius(d->atomRadius * d->hydrogenScale);
    }
    else{
        item->setRadius(d->atomRadius);
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
