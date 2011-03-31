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

#include "mmffpartialchargepredictor.h"

#include <chemkit/molecule.h>
#include <chemkit/pluginmanager.h>

#include "mmffplugin.h"
#include "mmffparameters.h"

// --- Construction and Destruction ---------------------------------------- //
MmffPartialChargePredictor::MmffPartialChargePredictor()
    : chemkit::PartialChargePredictor("mmff")
{
    m_typer = 0;
    m_parameters = 0;

    // load parameters
    const chemkit::Plugin *mmffPlugin = chemkit::PluginManager::instance()->plugin("mmff");
    if(!mmffPlugin){
        return;
    }

    m_parameters = new MmffParameters;
    bool ok = m_parameters->read(QString::fromStdString(mmffPlugin->dataPath()) + "mmff94.prm");
    if(!ok){
        delete m_parameters;
        m_parameters = 0;
    }
}

MmffPartialChargePredictor::~MmffPartialChargePredictor()
{
    delete m_parameters;
}

// --- Properties ---------------------------------------------------------- //
void MmffPartialChargePredictor::setAtomTyper(const MmffAtomTyper *typer)
{
    m_typer = typer;
}

// --- Partial Charges ----------------------------------------------------- //
chemkit::Float MmffPartialChargePredictor::partialCharge(int index) const
{
    return m_partialCharges.value(index, 0);
}

chemkit::Float MmffPartialChargePredictor::partialCharge(const chemkit::Atom *atom) const
{
    return partialCharge(atom->index());
}

void MmffPartialChargePredictor::assignPartialCharges(const chemkit::Molecule *molecule)
{
    if(!molecule){
        m_partialCharges.resize(0);
        return;
    }

    // ensure parameters are loaded
    if(!m_parameters){
        return;
    }

    // load atom typer
    const MmffAtomTyper *typer;
    if(!m_typer){
        typer = new MmffAtomTyper(molecule);
    }
    else{
        typer = m_typer;
    }

    // setup space for partial charges
    m_partialCharges.resize(molecule->size());
    m_partialCharges.fill(0);

    // assign partial charges for each atom
    for(int i = 0; i < molecule->size(); i++){
        const chemkit::Atom *atom = molecule->atom(i);
        int atomType = typer->typeNumber(atom);

        const MmffAtomParameters *atomParameters = m_parameters->atomParameters(atomType);
        if(!atomParameters){
            continue;
        }

        chemkit::Float q0 = typer->formalCharge(atom);
        chemkit::Float M = atomParameters->crd;
        chemkit::Float V = m_parameters->partialChargeParameters(atomType)->fcadj;
        chemkit::Float formalChargeSum = 0;
        chemkit::Float partialChargeSum = 0;

        if(V == 0){
            foreach(const chemkit::Atom *neighbor, atom->neighbors()){
                chemkit::Float neighborFormalCharge = typer->formalCharge(neighbor);

                if(neighborFormalCharge < 0){
                    q0 += neighborFormalCharge / (2 * neighbor->neighborCount());
                }
            }
        }

        if(atomType == 62){
            foreach(const chemkit::Atom *neighbor, atom->neighbors()){
                chemkit::Float neighborFormalCharge = typer->formalCharge(neighbor);

                if(neighborFormalCharge > 0){
                    q0 -= neighborFormalCharge / 2;
                }
            }
        }

        foreach(const chemkit::Atom *neighbor, atom->neighbors()){
            int neighborType = typer->typeNumber(neighbor);

            const MmffChargeParameters *chargeParameters = m_parameters->chargeParameters(atom, atomType, neighbor, neighborType);
            if(chargeParameters){
                partialChargeSum += -chargeParameters->bci;
            }
            else{
                chargeParameters = m_parameters->chargeParameters(neighbor, neighborType, atom, atomType);
                if(chargeParameters){
                    partialChargeSum += chargeParameters->bci;
                }
                else{
                    const MmffPartialChargeParameters *partialChargeParameters = m_parameters->partialChargeParameters(atomType);
                    const MmffPartialChargeParameters *neighborPartialChargeParameters = m_parameters->partialChargeParameters(neighborType);
                    if(!partialChargeParameters || !neighborPartialChargeParameters){
                        continue;
                    }

                    partialChargeSum += (partialChargeParameters->pbci - neighborPartialChargeParameters->pbci);
                }
            }

            formalChargeSum += typer->formalCharge(neighbor);
        }

        // equation 15 (p. 662)
        m_partialCharges[i] = (1 - M * V) * q0 + V * formalChargeSum + partialChargeSum;
    }

    // cleanup typer object
    if(!m_typer){
        delete typer;
    }
}
