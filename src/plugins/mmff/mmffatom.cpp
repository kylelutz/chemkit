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

#include "mmffatom.h"

#include "mmffforcefield.h"
#include "mmffparameters.h"

// --- Construction and Destruction ---------------------------------------- //
MmffAtom::MmffAtom(chemkit::ForceField *forceField, const chemkit::Atom *atom)
    : chemkit::ForceFieldAtom(forceField, atom)
{
    m_typeNumber = 0;
    m_formalCharge = 0;
}

// --- Properties ---------------------------------------------------------- //
const MmffForceField* MmffAtom::forceField() const
{
    return static_cast<const MmffForceField *>(chemkit::ForceFieldAtom::forceField());
}

void MmffAtom::setType(int typeNumber, chemkit::Float formalCharge)
{
    m_typeNumber = typeNumber;
    m_formalCharge = formalCharge;
}

QString MmffAtom::type() const
{
    return QString::number(m_typeNumber);
}

int MmffAtom::typeNumber() const
{
    return m_typeNumber;
}

chemkit::Float MmffAtom::formalCharge() const
{
    return m_formalCharge;
}

int MmffAtom::period() const
{
    return atom()->element().period();
}

bool MmffAtom::setCharge()
{
    const MmffParameters *parameters = forceField()->parameters();
    const MmffAtomParameters *atomParameters = this->parameters();
    if(!parameters || !atomParameters){
        return false;
    }

    chemkit::Float q0 = m_formalCharge;
    chemkit::Float M = atomParameters->crd;
    chemkit::Float V = parameters->partialChargeParameters(this)->fcadj;
    chemkit::Float formalChargeSum = 0;
    chemkit::Float partialChargeSum = 0;

    if(V == 0){
        foreach(const chemkit::Atom *neighborAtom, atom()->neighbors()){
            const MmffAtom *neighbor = forceField()->atom(neighborAtom);

            if(neighbor->formalCharge() < 0){
                q0 += neighbor->formalCharge() / (2 * neighborAtom->neighborCount());
            }
        }
    }

    if(typeNumber() == 62){
        foreach(const chemkit::Atom *neighborAtom, atom()->neighbors()){
            const MmffAtom *neighbor = forceField()->atom(neighborAtom);

            if(neighbor->formalCharge() > 0){
                q0 -= neighbor->formalCharge() / 2;
            }
        }
    }

    foreach(const chemkit::Atom *neighborAtom, atom()->neighbors()){
        const MmffAtom *neighbor = forceField()->atom(neighborAtom);

        const MmffChargeParameters *chargeParameters = parameters->chargeParameters(this, neighbor);
        if(chargeParameters){
            partialChargeSum += -chargeParameters->bci;
        }
        else{
            chargeParameters = parameters->chargeParameters(neighbor, this);
            if(chargeParameters){
                partialChargeSum += chargeParameters->bci;
            }
            else{
                const MmffPartialChargeParameters *partialChargeParameters = parameters->partialChargeParameters(this);
                const MmffPartialChargeParameters *neighborPartialChargeParameters = parameters->partialChargeParameters(neighbor);
                if(!partialChargeParameters || !neighborPartialChargeParameters){
                    return false;
                }

                partialChargeSum += (partialChargeParameters->pbci - neighborPartialChargeParameters->pbci);
            }
        }

        formalChargeSum += neighbor->formalCharge();
    }

    // equation 15 (p. 662)
    chemkit::Float charge = (1 - M * V) * q0 + V * formalChargeSum + partialChargeSum;

    ForceFieldAtom::setCharge(charge);
    return true;
}

// --- Parameters ---------------------------------------------------------- //
const MmffAtomParameters* MmffAtom::parameters() const
{
    const MmffParameters *mmffParameters = static_cast<const MmffForceField *>(forceField())->parameters();
    if(!mmffParameters)
        return 0;

    return mmffParameters->atomParameters(this);
}
