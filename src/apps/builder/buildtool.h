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
    void setAtomPosition(chemkit::Atom *atom, const chemkit::Point3 &position);
    chemkit::Bond* addBond(chemkit::Atom *a, chemkit::Atom *b, int order = 1);
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
