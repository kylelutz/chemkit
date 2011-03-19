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

#ifndef CHEMKIT_GRAPHICSMOLECULEITEM_H
#define CHEMKIT_GRAPHICSMOLECULEITEM_H

#include "graphics.h"

#include <chemkit/molecule.h>
#include <chemkit/moleculewatcher.h>

#include "graphicsitem.h"
#include "graphicsatomitem.h"
#include "graphicsbonditem.h"

namespace chemkit {

class GraphicsMoleculeItemPrivate;

class CHEMKIT_GRAPHICS_EXPORT GraphicsMoleculeItem : public QObject, public GraphicsItem
{
    Q_OBJECT

    public:
        // enumerations
        enum DisplayType {
            BallAndStick,
            Stick,
            SpaceFilling
        };

        // construction and destruction
        GraphicsMoleculeItem(const Molecule *molecule = 0);
        ~GraphicsMoleculeItem();

        // properties
        void setMolecule(const Molecule *molecule);
        const Molecule* molecule() const;
        void setDisplayType(DisplayType type);
        DisplayType displayType() const;
        void setAtomRadius(float radius);
        float atomRadius() const;
        void setBondRadius(float radius);
        float bondRadius() const;
        void setAtomVisible(const Atom *atom, bool visible);
        bool atomVisible(const Atom *atom) const;
        void setHydrogensVisible(bool visible);
        bool hydrogensVisible() const;
        void setHydrogenScale(float scale);
        float hydrogenScale() const;
        void setBondOrderVisible(bool visible);
        bool bondOrderVisible() const;
        void setAtomColoredBonds(bool atomColoredBonds);
        bool atomColoredBonds() const;

        // items
        GraphicsAtomItem* atomItem(const Atom *atom) const;
        GraphicsBondItem* bondItem(const Bond *bond) const;

        // painting
        virtual void paint(GraphicsPainter *painter);

        // static methods
        static QColor atomColor(int atomicNumber);
        static QColor atomColor(const Atom *atom);

    protected:
        // events
        void itemChanged(ItemChange change);

    private slots:
        // slots
        void atomAdded(const chemkit::Atom *atom);
        void atomRemoved(const chemkit::Atom *atom);
        void atomAtomicNumberChanged(const chemkit::Atom *atom);
        void atomPositionChanged(const chemkit::Atom *atom);
        void bondAdded(const chemkit::Bond *bond);
        void bondRemoved(const chemkit::Bond *bond);
        void bondOrderChanged(const chemkit::Bond *bond);

    private:
        GraphicsMoleculeItemPrivate* const d;
};

} // end chemkit namespace

#endif // CHEMKIT_GRAPHICSMOLECULEITEM_H
