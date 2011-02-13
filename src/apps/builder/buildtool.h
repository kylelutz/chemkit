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

#ifndef BUILDTOOL_H
#define BUILDTOOL_H

#include "buildertool.h"

class BuildTool : public QObject, public BuilderTool
{
    Q_OBJECT

    public:
        // construction and destruction
        BuildTool(BuilderWindow *builder);
        ~BuildTool();

        // properties
        void setElement(const chemkit::Element &element);
        chemkit::Element element() const;
        void setBondOrder(int bondOrder);
        int bondOrder() const;

        // settings
        virtual QWidget* settingsWidget();

        // events
        virtual void mousePressEvent(QMouseEvent *event);
        virtual void mouseReleaseEvent(QMouseEvent *event);
        virtual void mouseMoveEvent(QMouseEvent *event);


    private slots:
        void elementSelectorChanged(int index);
        void bondOrderSelectorChanged(int index);
        void addHydrogensChanged(int state);

    private:
        void beginMoleculeEdit();
        void endMoleculeEdit();
        chemkit::Atom* addAtom(int atomicNumber);
        void removeAtom(chemkit::Atom *atom);
        void setAtomAtomicNumber(chemkit::Atom *atom, int atomicNumber);
        void setAtomPosition(chemkit::Atom *atom, const chemkit::Point &position);
        chemkit::Bond* addBond(chemkit::Atom *a, chemkit::Atom *b, int order = chemkit::Bond::Single);
        void removeBond(chemkit::Bond *bond);
        void setBondOrder(chemkit::Bond *bond, int order);
        void adjustHydrogens(chemkit::Atom *atom);

    private:
        chemkit::Element m_element;
        int m_bondOrder;
        int m_intialElement;
        bool m_adjustHydrogens;
        QList<int> m_elements;
        QList<int> m_addedElements;
        chemkit::Atom *m_intialAtom;
        chemkit::Atom *m_movingAtom;
        chemkit::Atom *m_bondingAtom;
        chemkit::Bond *m_newBond;
        QComboBox *m_elementSelector;
        QComboBox *m_bondOrderSelector;
        QCheckBox *m_addHydrogensCheckBox;
        QSet<chemkit::Atom *> m_modifiedAtoms;
};

#endif // BUILDTOOL_H
