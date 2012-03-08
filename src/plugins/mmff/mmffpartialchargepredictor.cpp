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

#include "mmffpartialchargepredictor.h"

#include <chemkit/foreach.h>
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
    bool ok = m_parameters->read(mmffPlugin->dataPath() + "mmff94.prm");
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
void MmffPartialChargePredictor::setMolecule(const chemkit::Molecule *molecule)
{
    chemkit::PartialChargePredictor::setMolecule(molecule);

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
    std::fill(m_partialCharges.begin(), m_partialCharges.end(), 0);

    // assign partial charges for each atom
    for(size_t i = 0; i < molecule->size(); i++){
        const chemkit::Atom *atom = molecule->atom(i);
        int atomType = typer->typeNumber(atom);

        const MmffAtomParameters *atomParameters = m_parameters->atomParameters(atomType);
        if(!atomParameters){
            continue;
        }

        chemkit::Real q0 = typer->formalCharge(atom);
        chemkit::Real M = atomParameters->crd;
        chemkit::Real V = m_parameters->partialChargeParameters(atomType)->fcadj;
        chemkit::Real formalChargeSum = 0;
        chemkit::Real partialChargeSum = 0;

        if(V == 0){
            foreach(const chemkit::Atom *neighbor, atom->neighbors()){
                chemkit::Real neighborFormalCharge = typer->formalCharge(neighbor);

                if(neighborFormalCharge < 0){
                    q0 += neighborFormalCharge / (2 * neighbor->neighborCount());
                }
            }
        }

        if(atomType == 62){
            foreach(const chemkit::Atom *neighbor, atom->neighbors()){
                chemkit::Real neighborFormalCharge = typer->formalCharge(neighbor);

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

void MmffPartialChargePredictor::setAtomTyper(const MmffAtomTyper *typer)
{
    m_typer = typer;
}

// --- Partial Charges ----------------------------------------------------- //
chemkit::Real MmffPartialChargePredictor::partialCharge(const chemkit::Atom *atom) const
{
    return m_partialCharges[atom->index()];
}
