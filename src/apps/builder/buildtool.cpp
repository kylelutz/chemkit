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

#include "buildtool.h"

#include <chemkit/atom.h>
#include <chemkit/bond.h>
#include <chemkit/vector3.h>
#include <chemkit/element.h>
#include <chemkit/periodictabledialog.h>

// --- Construction and Destruction ---------------------------------------- //
BuildTool::BuildTool(BuilderWindow *builder)
    : QObject(),
      BuilderTool(builder)
{
    m_element = chemkit::Atom::Carbon;
    m_bondOrder = chemkit::Bond::Single;
    m_adjustHydrogens = true;

    m_intialAtom = 0;
    m_movingAtom = 0;
    m_bondingAtom = 0;
    m_newBond = 0;

    // elements to display in element selector
    m_elements.append(chemkit::Atom::Carbon);
    m_elements.append(chemkit::Atom::Nitrogen);
    m_elements.append(chemkit::Atom::Oxygen);
    m_elements.append(chemkit::Atom::Chlorine);
    m_elements.append(chemkit::Atom::Bromine);
    m_elements.append(chemkit::Atom::Hydrogen);
    m_elements.append(chemkit::Atom::Phosphorus);
    m_elements.append(chemkit::Atom::Sulfur);
}

BuildTool::~BuildTool()
{
}

/// Sets the current element to \p element.
void BuildTool::setElement(const chemkit::Element &element)
{
    if(element.isValid()){
        m_element = element;
    }

    if(m_elementSelector){
        if(m_elements.contains(element.atomicNumber())){
            m_elementSelector->setCurrentIndex(m_elements.indexOf(element.atomicNumber()));
        }
        else if(m_addedElements.contains(element.atomicNumber())){
            m_elementSelector->setCurrentIndex(m_elementSelector->findText(element.name().c_str()));
        }
        else{
            m_elementSelector->removeItem(m_elementSelector->count() - 1);
            m_elementSelector->addItem(element.name().c_str(), element.atomicNumber());
            m_elementSelector->addItem("Other", -1);
            m_elementSelector->update();
            m_addedElements.append(element.atomicNumber());

            m_elementSelector->setCurrentIndex(m_elementSelector->count() - 2);
        }
    }
}

/// Returns the current element.
chemkit::Element BuildTool::element() const
{
    return m_element;
}

/// Sets the current bond order to \p bondOrder.
void BuildTool::setBondOrder(int bondOrder)
{
    m_bondOrder = bondOrder;
}

/// Returns the current bond order.
int BuildTool::bondOrder() const
{
    return m_bondOrder;
}

// --- Settings ------------------------------------------------------------ //
QWidget* BuildTool::settingsWidget()
{
    QWidget *widget = new QWidget;

    QFormLayout *layout = new QFormLayout;

    // element selector
    m_elementSelector = new QComboBox;
    foreach(int element, m_elements){
        m_elementSelector->addItem(chemkit::Element(element).name().c_str(), element);
    }
    if(!m_addedElements.isEmpty()){
        foreach(int element, m_addedElements){
            m_elementSelector->addItem(chemkit::Element(element).name().c_str(), element);
        }
    }
    m_elementSelector->addItem("Other", -1);
    connect(m_elementSelector, SIGNAL(currentIndexChanged(int)), SLOT(elementSelectorChanged(int)));
    layout->addRow("Element:", m_elementSelector);

    // bond order selector
    m_bondOrderSelector = new QComboBox;
    m_bondOrderSelector->addItem("Single");
    m_bondOrderSelector->addItem("Double");
    m_bondOrderSelector->addItem("Triple");
    m_bondOrderSelector->setCurrentIndex(m_bondOrder - 1);
    connect(m_bondOrderSelector, SIGNAL(currentIndexChanged(int)), SLOT(bondOrderSelectorChanged(int)));
    layout->addRow("Bond Order:", m_bondOrderSelector);

    // add hydrogens checkbox
    m_addHydrogensCheckBox = new QCheckBox("Auto Add Hydrogens");
    m_addHydrogensCheckBox->setChecked(m_adjustHydrogens);
    connect(m_addHydrogensCheckBox, SIGNAL(stateChanged(int)), SLOT(addHydrogensChanged(int)));
    layout->addRow(m_addHydrogensCheckBox);

    widget->setLayout(layout);

    return widget;
}

// --- Events -------------------------------------------------------------- //
void BuildTool::mousePressEvent(QMouseEvent *event)
{
    beginMoleculeEdit();

    // left button
    if(event->button() == Qt::LeftButton){
        chemkit::GraphicsItem *item = view()->itemAt(event->x(), event->y());

        // item under cursor
        if(item){
            // atom under cursor
            if(item->type() == chemkit::GraphicsItem::AtomItem){
                chemkit::GraphicsAtomItem *atomItem = static_cast<chemkit::GraphicsAtomItem *>(item);
                chemkit::Atom *atom = const_cast<chemkit::Atom *>(atomItem->atom());
                m_intialElement = atom->atomicNumber();

                if(atom->atomicNumber() != m_element.atomicNumber()){
                    setAtomAtomicNumber(atom, m_element.atomicNumber());
                }

                m_intialAtom = atom;
            }

            // bond under cursor
            else if(item->type() == chemkit::GraphicsItem::BondItem){
                chemkit::GraphicsBondItem *bondItem = static_cast<chemkit::GraphicsBondItem *>(item);
                chemkit::Bond *bond = const_cast<chemkit::Bond *>(bondItem->bond());
                setBondOrder(bond, (bond->order() % 3) + 1);
            }
        }

        // no item under cursor
        else{
            // add new atom
            chemkit::Atom *atom = addAtom(m_element.atomicNumber());
            chemkit::Point3f position = view()->unproject(event->x(), event->y(), editor()->molecule()->center().cast<float>());
            setAtomPosition(atom, position.cast<chemkit::Real>());
            m_intialAtom = atom;
            m_intialElement = m_element.atomicNumber();
        }

        m_movingAtom = 0;
        m_bondingAtom = 0;
        m_newBond = 0;
    }

    // right button
    else if(event->button() == Qt::RightButton){
        chemkit::GraphicsItem *item = view()->itemAt(event->x(), event->y());

        if(item){
            if(item->type() == chemkit::GraphicsItem::AtomItem){
                chemkit::GraphicsAtomItem *atomItem = static_cast<chemkit::GraphicsAtomItem *>(item);
                chemkit::Atom *atom = const_cast<chemkit::Atom *>(atomItem->atom());
                removeAtom(atom);
            }
            else if(item->type() == chemkit::GraphicsItem::BondItem){
                chemkit::GraphicsBondItem *bondItem = static_cast<chemkit::GraphicsBondItem *>(item);
                chemkit::Bond *bond = const_cast<chemkit::Bond *>(bondItem->bond());
                removeBond(bond);
            }
        }
    }

    view()->update();
}

void BuildTool::mouseMoveEvent(QMouseEvent *event)
{
    if(!m_intialAtom){
        return;
    }

    chemkit::GraphicsItem *item = 0;

    if(m_movingAtom){
        std::vector<chemkit::GraphicsItem *> items = view()->itemsAt(event->x(), event->y());

        if(items.empty()){
            item = 0;
        }
        else if(items.size() == 1){
            item = items[0];
        }
        else{
            chemkit::GraphicsItem *nearestItem = items[0];

            // if the nearest item that we intersect is the atom item for
            // the currently moving atom, then select the next item beneath it.
            if(nearestItem->type() == chemkit::GraphicsItem::AtomItem &&
               static_cast<chemkit::GraphicsAtomItem *>(nearestItem)->atom() == m_movingAtom){
                item = items[1];
            }
            else{
                item = items[0];
            }
        }
    }
    else{
        item = view()->itemAt(event->x(), event->y());
    }

    // cursor over nothing
    if(!item){
        if(!m_movingAtom){
            setAtomAtomicNumber(m_intialAtom, m_intialElement);
            m_movingAtom = addAtom(m_element.atomicNumber());
            addBond(m_intialAtom, m_movingAtom, bondOrder());
            chemkit::Point3f position = view()->unproject(event->x(), event->y(), m_intialAtom->position().cast<float>());
            setAtomPosition(m_movingAtom, position.cast<chemkit::Real>());

            if(m_newBond){
                removeBond(m_newBond);
                m_newBond = 0;
                m_bondingAtom = 0;
            }
        }
        else{
            chemkit::Point3f newPosition = view()->unproject(event->x(), event->y(), m_movingAtom->position().cast<float>());
            setAtomPosition(m_movingAtom, newPosition.cast<chemkit::Real>());
        }
    }

    // cursor over atom item
    else if(item->type() == chemkit::GraphicsItem::AtomItem){
        chemkit::GraphicsAtomItem *atomItem = static_cast<chemkit::GraphicsAtomItem *>(item);
        chemkit::Atom *atom = const_cast<chemkit::Atom *>(atomItem->atom());

        // over initial atom
        if(atom == m_intialAtom){
            if(m_movingAtom){
                removeAtom(m_movingAtom);
                m_movingAtom = 0;
                setAtomAtomicNumber(m_intialAtom, m_element.atomicNumber());
            }
        }
        // over moving atom
        else if(atom == m_movingAtom){
            chemkit::Point3f newPosition = view()->unproject(event->x(), event->y(), m_movingAtom->position().cast<float>());
            setAtomPosition(m_movingAtom, newPosition.cast<chemkit::Real>());
        }
        // over new atom
        else{
            if(m_movingAtom){
                removeAtom(m_movingAtom);
                m_movingAtom = 0;
            }

            if(m_newBond && atom != m_bondingAtom){
                removeBond(m_newBond);
                m_newBond = 0;
                m_bondingAtom = 0;
            }

            if(!m_intialAtom->isBondedTo(atom)){
                m_newBond = addBond(atom, m_intialAtom);
                setBondOrder(m_newBond, m_bondOrder);
                m_bondingAtom = atom;
            }
        }
    }

    view()->update();
}

void BuildTool::mouseReleaseEvent(QMouseEvent *event)
{
    if(event->button() == Qt::LeftButton){
        m_intialAtom = 0;
        m_movingAtom = 0;
        m_bondingAtom = 0;
        m_newBond = 0;
    }

    if(m_adjustHydrogens){
        foreach(chemkit::Atom *atom, m_modifiedAtoms){
            adjustHydrogens(atom);
        }
    }
    m_modifiedAtoms.clear();

    endMoleculeEdit();
}

// --- Slots --------------------------------------------------------------- //
void BuildTool::elementSelectorChanged(int index)
{
    int atomicNumber = m_elementSelector->itemData(index).toInt();
    if(atomicNumber == -1){
        const chemkit::Element &element = chemkit::PeriodicTableDialog::getElement(builder(), "Select Element");

        if(element.isValid()){
            setElement(element);
        }
        else{
            setElement(chemkit::Atom::Carbon);
        }
    }
    else{
        setElement(atomicNumber);
    }
}

void BuildTool::bondOrderSelectorChanged(int index)
{
    setBondOrder(index + 1);
}

void BuildTool::addHydrogensChanged(int state)
{
    if(state == Qt::Checked){
        m_adjustHydrogens = true;
    }
    else{
        m_adjustHydrogens = false;
    }
}

// --- Internal Methods ---------------------------------------------------- //
void BuildTool::beginMoleculeEdit()
{
    builder()->beginMoleculeEdit();
}

void BuildTool::endMoleculeEdit()
{
    // do hydrogen adjustment
    if(m_adjustHydrogens){
        foreach(chemkit::Atom *atom, m_modifiedAtoms){
            adjustHydrogens(atom);
        }
    }
    m_modifiedAtoms.clear();

    builder()->endMoleculeEdit();
}

chemkit::Atom* BuildTool::addAtom(int atomicNumber)
{
    chemkit::Atom *atom = editor()->addAtom(atomicNumber);
    m_modifiedAtoms.insert(atom);
    return atom;
}

void BuildTool::removeAtom(chemkit::Atom *atom)
{
    foreach(chemkit::Atom *neighbor, atom->neighbors()){
        m_modifiedAtoms.insert(neighbor);
    }

    editor()->removeAtom(atom);
    m_modifiedAtoms.remove(atom);
}

void BuildTool::setAtomAtomicNumber(chemkit::Atom *atom, int atomicNumber)
{
    m_modifiedAtoms.insert(atom);
    editor()->setAtomElement(atom, atomicNumber);
}

void BuildTool::setAtomPosition(chemkit::Atom *atom, const chemkit::Point3 &position)
{
    editor()->setAtomPosition(atom, position);
}

chemkit::Bond* BuildTool::addBond(chemkit::Atom *a, chemkit::Atom *b, int order)
{
    m_modifiedAtoms.insert(a);
    m_modifiedAtoms.insert(b);
    return editor()->addBond(a, b, order);
}

void BuildTool::removeBond(chemkit::Bond *bond)
{
    m_modifiedAtoms.insert(bond->atom1());
    m_modifiedAtoms.insert(bond->atom2());
    editor()->removeBond(bond);
}

void BuildTool::setBondOrder(chemkit::Bond *bond, int order)
{
    m_modifiedAtoms.insert(bond->atom1());
    m_modifiedAtoms.insert(bond->atom2());
    editor()->setBondOrder(bond, order);
}

void BuildTool::adjustHydrogens(chemkit::Atom *atom)
{
    // remove lone hydrogens
    if(atom->is(chemkit::Atom::Hydrogen) && atom->neighborCount() < 2){
        editor()->removeAtom(atom);
        m_modifiedAtoms.remove(atom);
        return;
    }

    // add hydrogens
    while(atom->valence() < atom->expectedValence()){
        chemkit::Atom *hydrogen = editor()->addAtom(chemkit::Atom::Hydrogen);
        editor()->setAtomPosition(hydrogen, atom->position() + chemkit::Vector3::Random().normalized());
        editor()->addBond(atom, hydrogen);
    }

    // remove hydrogens
    while(atom->valence() > atom->expectedValence()){
        bool foundTerminalHydrogen = false;

        foreach(chemkit::Atom *neighbor, atom->neighbors()){
            if(neighbor->isTerminalHydrogen()){
                editor()->removeAtom(neighbor);
                m_modifiedAtoms.remove(neighbor);
                foundTerminalHydrogen = true;
                break;
            }
        }

        if(!foundTerminalHydrogen){
            // no more hydrogens to remove
            break;
        }
    }
}
