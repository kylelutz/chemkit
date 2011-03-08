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

#include "gasteigerpartialchargepredictor.h"

#include <chemkit/atom.h>
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
GasteigerPartialChargePredictor::GasteigerPartialChargePredictor()
    : chemkit::PartialChargePredictor("gasteiger")
{
}

GasteigerPartialChargePredictor::~GasteigerPartialChargePredictor()
{
}

// --- Partial Charges ----------------------------------------------------- //
chemkit::Float GasteigerPartialChargePredictor::partialCharge(int index) const
{
    return m_charges.value(index, 0);
}

void GasteigerPartialChargePredictor::assignPartialCharges(const chemkit::Molecule *molecule)
{
    if(!molecule){
        m_charges.resize(0);
        return;
    }

    m_charges.resize(molecule->atomCount());
    m_electronegativies.resize(molecule->atomCount());
    m_parameters.resize(molecule->atomCount());

    // initialize charges, electronegativities and parameters
    for(int i = 0; i < molecule->atomCount(); i++){
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
        for(int i = 0; i < molecule->atomCount(); i++){
            const chemkit::Atom *atom = molecule->atom(i);

            chemkit::Float qi = 0;
            chemkit::Float Xi = m_electronegativies[i];
            const GasteigerParameters *pi = m_parameters[i];

            foreach(const chemkit::Atom *neighbor, atom->neighbors()){
                int j = neighbor->index();
                chemkit::Float Xj = m_electronegativies[j];
                const GasteigerParameters *pj = m_parameters[j];

                chemkit::Float scale;
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
        for(int i = 0; i < molecule->atomCount(); i++){
            const GasteigerParameters *pi = m_parameters[i];
            chemkit::Float Qi = m_charges[i];

            m_electronegativies[i] = pi->a + pi->b * Qi + pi->c * pow(Qi, 2);
        }
    }
}

// --- Internal Methods ---------------------------------------------------- //
const GasteigerParameters* GasteigerPartialChargePredictor::atomParameters(const chemkit::Atom *atom) const
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
