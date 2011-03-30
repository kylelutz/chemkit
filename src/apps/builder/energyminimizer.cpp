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

#include "energyminimizer.h"

#include <QtConcurrentRun>

// --- Construction and Destruction ---------------------------------------- //
EnergyMinimizer::EnergyMinimizer(chemkit::Molecule *molecule)
    : QObject()
{
    m_molecule = molecule;
    m_moleculeChanged = true;
    m_state = Stopped;
    m_forceField = 0;
    m_forceFieldName = "uff";
    connect(&m_minimizationWatcher, SIGNAL(finished()), SLOT(minimizationStepFinished()));
}

EnergyMinimizer::~EnergyMinimizer()
{
    delete m_forceField;
}

// --- Properties ---------------------------------------------------------- //
void EnergyMinimizer::setMolecule(chemkit::Molecule *molecule)
{
    if(molecule == m_molecule){
        return;
    }

    m_molecule = molecule;
    m_moleculeChanged = true;
}

chemkit::Molecule* EnergyMinimizer::molecule() const
{
    return m_molecule;
}

void EnergyMinimizer::setMoleculeChanged(bool changed)
{
    m_moleculeChanged = changed;
}

bool EnergyMinimizer::moleculeChanged() const
{
    return m_moleculeChanged;
}

void EnergyMinimizer::setForceField(const QString &name)
{
    m_forceFieldName = name;
    m_moleculeChanged = true;
}

chemkit::ForceField* EnergyMinimizer::forceField() const
{
    return m_forceField;
}

int EnergyMinimizer::state() const
{
    return m_state;
}

QString EnergyMinimizer::stateString() const
{
    switch(m_state){
        case Running: return "Running";
        case Stopped: return "Stopped";
        case SettingUp: return "Setting Up";
        case SetupFailed: return "Setup Failed";
        case UpdateReady: return "Update Ready";
        case Converged: return "Converged";
        default: return "Unknown";
    }
}

// --- Optimization -------------------------------------------------------- //
chemkit::Float EnergyMinimizer::energy() const
{
    if(m_forceField)
        return m_forceField->energy();
    else
        return 0;
}

// --- Slots --------------------------------------------------------------- //
void EnergyMinimizer::start()
{
    if(!m_molecule || m_molecule->isEmpty()){
        setState(SetupFailed);
        return;
    }

    if(m_moleculeChanged){
        delete m_forceField;

        m_forceField = chemkit::ForceField::create(m_forceFieldName.toStdString());
        if(!m_forceField){
            setState(SetupFailed);
            return;
        }

        setState(SettingUp);

        m_forceField->addMolecule(m_molecule);
        m_forceField->setup();

        if(!m_forceField->isSetup()){
            setState(SetupFailed);
            return;
        }

        m_moleculeChanged = false;
    }

    QFuture<bool> future = m_forceField->minimizationStepAsync();
    m_minimizationWatcher.setFuture(future);

    setState(Running);
}

void EnergyMinimizer::stop()
{
    setState(Stopped);
}

// --- Internal Methods ---------------------------------------------------- //
void EnergyMinimizer::setState(int state)
{
    if(state == m_state)
        return;

    m_state = state;
    emit stateChanged(state);
}

void EnergyMinimizer::minimizationStepFinished()
{
    bool converged = m_minimizationWatcher.result();

    if(m_state == Stopped){
        return;
    }

    if(converged){
        setState(Converged);
    }
    else{
        setState(UpdateReady);
    }
}
