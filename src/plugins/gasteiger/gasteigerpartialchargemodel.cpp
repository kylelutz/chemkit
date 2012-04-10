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

#include "gasteigerpartialchargemodel.h"

#include <chemkit/atom.h>
#include <chemkit/foreach.h>
#include <chemkit/molecule.h>

namespace {

// Parameters from Table 1 in [Gasteiger 1980]
const GasteigerParameters Parameters[] = {
    {7.17, 6.24, -0.56}, // 0 - H
    {7.98, 9.18, 1.88}, // 1 - C (sp3)
    {8.79, 9.32, 1.51}, // 2 - C (sp2)
    {10.39, 9.45, 0.73}, // 3 - C (sp)
    {11.54, 10.82, 1.36}, // 4 - N (sp3)
    {12.87, 11.15, 0.85}, // 5 - N (sp2)
    {15.68, 11.7, -0.27}, // 6 - N (sp)
    {14.18, 12.92, 1.39}, // 7 - O (sp3)
    {17.07, 13.79, 0.47}, // 8 - O (sp2)
    {14.66, 13.85, 2.31}, // 9 - F
    {11.00, 9.69, 1.35}, // 10 - Cl
    {10.08, 8.47, 1.16}, // 11 - Br
    {9.90, 7.96, 0.96}, // 12 - I
    {10.14, 9.13, 1.38} // 13 - S
};

} // end anonymous namespace

// --- Construction and Destruction ---------------------------------------- //
GasteigerPartialChargeModel::GasteigerPartialChargeModel()
    : chemkit::PartialChargeModel("gasteiger")
{
}

// --- Properties ---------------------------------------------------------- //
void GasteigerPartialChargeModel::setMolecule(const chemkit::Molecule *molecule)
{
    chemkit::PartialChargeModel::setMolecule(molecule);

    if(!molecule){
        m_charges.resize(0);
        return;
    }

    m_charges.resize(molecule->atomCount());
    m_electronegativies.resize(molecule->atomCount());
    m_parameters.resize(molecule->atomCount());

    // initialize charges, electronegativities and parameters
    for(size_t i = 0; i < molecule->atomCount(); i++){
        const chemkit::Atom *atom = molecule->atom(i);

        const GasteigerParameters *parameters = atomParameters(atom);
        if(!parameters){
            return;
        }

        m_parameters[i] = parameters;
        m_charges[i] = 0;
        m_electronegativies[i] = parameters->a;
    }

    // run algorithm for six iterations
    for(int iteration = 1; iteration <= 6; iteration++){

        // calculate charges
        for(size_t i = 0; i < molecule->atomCount(); i++){
            const chemkit::Atom *atom = molecule->atom(i);

            chemkit::Real qi = 0;
            chemkit::Real Xi = m_electronegativies[i];
            const GasteigerParameters *pi = m_parameters[i];

            foreach(const chemkit::Atom *neighbor, atom->neighbors()){
                int j = neighbor->index();
                chemkit::Real Xj = m_electronegativies[j];
                const GasteigerParameters *pj = m_parameters[j];

                chemkit::Real scale;
                if(Xj > Xi){
                    if(atom->is(chemkit::Atom::Hydrogen)){
                        scale = 1.0 / 20.02;
                    }
                    else{
                        scale = 1.0 / (pi->a + pi->b + pi->c);
                    }
                }
                else{
                    if(neighbor->is(chemkit::Atom::Hydrogen)){
                        scale = 1.0 / 20.02;
                    }
                    else{
                        scale = 1.0 / (pj->a + pj->b + pj->c);
                    }
                }

                qi += scale * (Xj - Xi);
            }

            m_charges[i] += qi * pow(0.5, iteration);
        }

        // calculate electronegativities
        for(size_t i = 0; i < molecule->atomCount(); i++){
            const GasteigerParameters *pi = m_parameters[i];
            chemkit::Real Qi = m_charges[i];

            m_electronegativies[i] = pi->a + pi->b * Qi + pi->c * pow(Qi, 2);
        }
    }
}

// --- Partial Charges ----------------------------------------------------- //
chemkit::Real GasteigerPartialChargeModel::partialCharge(const chemkit::Atom *atom) const
{
    return m_charges[atom->index()];
}

// --- Internal Methods ---------------------------------------------------- //
const GasteigerParameters* GasteigerPartialChargeModel::atomParameters(const chemkit::Atom *atom) const
{
    if(atom->is(chemkit::Atom::Hydrogen)){
        return &Parameters[0];
    }
    else if(atom->is(chemkit::Atom::Carbon)){
        if(atom->neighborCount() == 4){
            return &Parameters[1];
        }
        else if(atom->neighborCount() == 3){
            return &Parameters[2];
        }
        else if(atom->neighborCount() == 2){
            return &Parameters[3];
        }
    }
    else if(atom->is(chemkit::Atom::Nitrogen)){
        if(atom->neighborCount() == 3){
            return &Parameters[4];
        }
        else if(atom->neighborCount() == 2){
            return &Parameters[5];
        }
        else if(atom->neighborCount() == 1){
            return &Parameters[6];
        }
    }
    else if(atom->is(chemkit::Atom::Oxygen)){
        if(atom->neighborCount() == 2){
            return &Parameters[7];
        }
        else if(atom->neighborCount() == 1){
            return &Parameters[8];
        }
    }
    else if(atom->is(chemkit::Atom::Fluorine)){
        return &Parameters[9];
    }
    else if(atom->is(chemkit::Atom::Chlorine)){
        return &Parameters[10];
    }
    else if(atom->is(chemkit::Atom::Bromine)){
        return &Parameters[11];
    }
    else if(atom->is(chemkit::Atom::Iodine)){
        return &Parameters[12];
    }
    else if(atom->is(chemkit::Atom::Sulfur)){
        return &Parameters[13];
    }

    return 0;
}
