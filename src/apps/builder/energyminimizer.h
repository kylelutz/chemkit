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

#ifndef GEOMETRYOPTIMIZER_H
#define GEOMETRYOPTIMIZER_H

#include <QObject>
#include <QFuture>

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
            Converged,
        };

        // construction and destruction
        EnergyMinimizer(chemkit::Molecule *molecule = 0);
        ~EnergyMinimizer();

        // properties
        void setMolecule(chemkit::Molecule *molecule);
        chemkit::Molecule* molecule() const;
        void setForceField(const QString &name);
        chemkit::ForceField* forceField() const;
        int state() const;
        QString stateString() const;

        // optimization
        chemkit::Float energy() const;
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
        chemkit::Molecule *m_molecule;
        chemkit::ForceField *m_forceField;
        QString m_forceFieldName;
        int m_state;
        QFutureWatcher<bool> m_minimizationWatcher;
};

#endif // GEOMETRYOPTIMIZER_H
