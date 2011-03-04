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

// This file implements the OPLS force field. See [Jorgensen 1996].

#include "oplsforcefield.h"

#include <chemkit/plugin.h>
#include <chemkit/molecule.h>
#include <chemkit/pluginmanager.h>
#include <chemkit/forcefieldinteractions.h>

#include "oplsatomtyper.h"
#include "oplsparameters.h"
#include "oplscalculation.h"

// --- Construction and Destruction ---------------------------------------- //
OplsForceField::OplsForceField()
    : chemkit::ForceField("opls")
{
    m_parameters = 0;
    setFlags(chemkit::ForceField::AnalyticalGradient);

    const chemkit::Plugin *oplsPlugin = chemkit::PluginManager::instance()->plugin("opls");
    if(oplsPlugin){
        m_parameters = new OplsParameters(oplsPlugin->dataPath() + "oplsaa.prm");
    }
}

OplsForceField::~OplsForceField()
{
    delete m_parameters;
}

// --- Parameterization ---------------------------------------------------- //
bool OplsForceField::setup()
{
    bool failed = false;

    foreach(const chemkit::Molecule *molecule, molecules()){
        if(!setupMolecule(molecule)){
            failed = true;
        }
    }

    if(!m_parameters){
        return false;
    }

    foreach(chemkit::ForceFieldCalculation *calculation, calculations()){
        bool setup = static_cast<OplsCalculation *>(calculation)->setup(m_parameters);

        if(!setup){
            failed = true;
        }

        setCalculationSetup(calculation, setup);
    }

    return !failed;
}

bool OplsForceField::setupMolecule(const chemkit::Molecule *molecule)
{
    OplsAtomTyper typer(molecule);

    foreach(const chemkit::Atom *atom, molecule->atoms()){
        chemkit::ForceFieldAtom *forceFieldAtom = new chemkit::ForceFieldAtom(this, atom);
        forceFieldAtom->setType(typer.typeString(atom));
        addAtom(forceFieldAtom);
    }

    chemkit::ForceFieldInteractions interactions(molecule, this);

    // bond strech pairs
    QPair<const chemkit::ForceFieldAtom *, const chemkit::ForceFieldAtom *> bondedPair;
    foreach(bondedPair, interactions.bondedPairs()){
        addCalculation(new OplsBondStrechCalculation(bondedPair.first, bondedPair.second));
    }

    // angle bend groups
    QVector<const chemkit::ForceFieldAtom *> angleGroup;
    foreach(angleGroup, interactions.angleGroups()){
        addCalculation(new OplsAngleBendCalculation(angleGroup[0], angleGroup[1], angleGroup[2]));
    }

    // torsion groups
    QVector<const chemkit::ForceFieldAtom *> torsionGroup;
    foreach(torsionGroup, interactions.torsionGroups()){
        addCalculation(new OplsTorsionCalculation(torsionGroup[0], torsionGroup[1], torsionGroup[2], torsionGroup[3]));
    }

    // nonbonded pairs
    QPair<const chemkit::ForceFieldAtom *, const chemkit::ForceFieldAtom *> nonbondedPair;
    foreach(nonbondedPair, interactions.nonbondedPairs()){
        addCalculation(new OplsNonbondedCalculation(nonbondedPair.first, nonbondedPair.second));
    }

    return true;
}
