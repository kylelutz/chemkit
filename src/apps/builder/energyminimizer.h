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

#ifndef GEOMETRYOPTIMIZER_H
#define GEOMETRYOPTIMIZER_H

#include <QObject>
#include <QFuture>
#include <QFutureWatcher>

#include <chemkit/molecule.h>
#include <chemkit/forcefield.h>

class EnergyMinimizer : public QObject
{
    Q_OBJECT

public:
    // enumerations
    enum State {
        Running,
        Stopped,
        SettingUp,
        SetupFailed,
        UpdateReady,
        Converged
    };

    // construction and destruction
    EnergyMinimizer(chemkit::Molecule *molecule = 0);
    ~EnergyMinimizer();

    // properties
    void setMolecule(chemkit::Molecule *molecule);
    chemkit::Molecule* molecule() const;
    void setMoleculeChanged(bool changed);
    bool moleculeChanged() const;
    void setForceField(const QString &name);
    chemkit::ForceField* forceField() const;
    int state() const;
    QString stateString() const;

    // optimization
    chemkit::Real energy() const;
    void reload() const;

public slots:
    void start();
    void stop();

signals:
    void forceFieldChanged(const chemkit::ForceField *forceField);
    void stateChanged(int state);

private slots:
    void minimizationStepFinished();

private:
    void setState(int state);

private:
    bool m_moleculeChanged;
    chemkit::Molecule *m_molecule;
    chemkit::ForceField *m_forceField;
    QString m_forceFieldName;
    int m_state;
    QFutureWatcher<bool> m_minimizationWatcher;
};

#endif // GEOMETRYOPTIMIZER_H
