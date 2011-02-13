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

#include <chemkit/molecule.h>

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

void MmffAtom::setType()
{
    const chemkit::Atom *atom = this->atom();

    switch(atom->atomicNumber()){
        // carbon
        case chemkit::Atom::Carbon:
            setCarbonType();
            break;

        // nitrogen
        case chemkit::Atom::Nitrogen:
            setNitrogenType();
            break;

        // oxygen
        case chemkit::Atom::Oxygen:
            setOxygenType();
            break;

        // phosphorus
        case chemkit::Atom::Phosphorus:
            if(atom->neighborCount() == 4){
                setType(25);
            }
            else if(atom->neighborCount() == 3){
                if(atom->isBondedTo(chemkit::Atom::Carbon, chemkit::Bond::Double))
                    setType(75);
                else
                    setType(26);
            }
            else if(atom->neighborCount() == 2 && atom->isBondedTo(chemkit::Atom::Carbon)){
                setType(75);
            }
            break;

        // sulfur
        case chemkit::Atom::Sulfur:
            setSulfurType();
            break;

        // fluorine
        case chemkit::Atom::Fluorine:
            if(atom->valence() > 0){
                setType(11);
            }
            else{
                setType(89, -1.0);
            }
            break;

        // chlorine
        case chemkit::Atom::Chlorine:
            if(atom->neighborCount(chemkit::Atom::Oxygen) == 4){
                setType(77);
            }
            else if(atom->valence() > 0){
                setType(12);
            }
            else{
                setType(90, -1.0);
            }
            break;

        // bromine
        case chemkit::Atom::Bromine:
            if(atom->valence() > 0){
                setType(13);
            }
            else{
                setType(91, -1.0);
            }
            break;

        // iodine
        case chemkit::Atom::Iodine:
            if(atom->valence() > 0){
                setType(14);
            }
            break;

        // iron
        case chemkit::Atom::Iron:
            if(qRound(atom->partialCharge()) == 2){
                setType(87, 2.0);
            }
            else{
                setType(88, 3.0);
            }
            break;

        // lithium
        case chemkit::Atom::Lithium:
            setType(92, 1.0);
            break;

        // sodium
        case chemkit::Atom::Sodium:
            setType(93, 1.0);
            break;

        // potassium
        case chemkit::Atom::Potassium:
            setType(94, 1.0);
            break;

        // zinc
        case chemkit::Atom::Zinc:
            setType(95, 2.0);
            break;

        // calcium
        case chemkit::Atom::Calcium:
            setType(96, 2.0);
            break;

        // copper
        case chemkit::Atom::Copper:
            if(qRound(atom->partialCharge()) == 2){
                setType(98, 2.0);
            }
            else{
                setType(97, 1.0);
            }
            break;

        // magnesium
        case chemkit::Atom::Magnesium:
            setType(99, 2.0);
            break;

        // silicon
        case chemkit::Atom::Silicon:
            setType(19);
            break;

        default:
            break;
    }
}

void MmffAtom::setHydrogenType(const MmffAtom *neighborAtom)
{
    Q_ASSERT(neighborAtom);
    Q_ASSERT(atom()->isTerminalHydrogen());

    const chemkit::Atom *atom = this->atom();
    const chemkit::Atom *neighbor = neighborAtom->atom();
    int neighborType = neighborAtom->typeNumber();

    // carbon
    if(neighbor->is(chemkit::Atom::Carbon)){
        setType(5);
    }

    // nitrogen
    else if(neighbor->is(chemkit::Atom::Nitrogen)){
        if(neighborType == 8 || neighborType == 39 || neighborType == 45 ||
           neighborType == 62 || neighborType == 67 || neighborType == 68){
            setType(23);
        }
        else if(neighborType == 9){
            setType(27);
        }
        else if(neighborType == 10 || neighborType == 40 || neighborType == 42 ||
                neighborType == 43 || neighborType == 48){
            setType(28);
        }
        else if(neighborType == 34 || neighborType == 54 || neighborType == 55 ||
                neighborType == 56 || neighborType == 58 || neighborType == 81){
            setType(36);
        }
    }

    // oxygen
    else if(neighbor->is(chemkit::Atom::Oxygen)){
        if(neighbor->isBondedTo(chemkit::Atom::Sulfur)){
            setType(33);
        }
        else if(neighborType == 6){
            bool imineOrEnol = false;
            bool carboxylicAcid = false;
            bool phosphate = false;

            foreach(const chemkit::Atom *secondNeighbor, neighbor->neighbors()){
                if(secondNeighbor == atom)
                    continue;

                if((secondNeighbor->is(chemkit::Atom::Carbon) || secondNeighbor->is(chemkit::Atom::Phosphorus)) &&
                   secondNeighbor->isBondedTo(chemkit::Atom::Oxygen, chemkit::Bond::Double)){
                    carboxylicAcid = true;
                    break;
                }
                else if(secondNeighbor->is(chemkit::Atom::Carbon) &&
                        (secondNeighbor->isBondedTo(chemkit::Atom::Carbon, chemkit::Bond::Double) ||
                        secondNeighbor->isBondedTo(chemkit::Atom::Nitrogen, chemkit::Bond::Double))){
                    imineOrEnol = true;
                    break;
                }
                else if(secondNeighbor->is(chemkit::Atom::Phosphorus) &&
                        secondNeighbor->neighborCount(chemkit::Atom::Oxygen) >= 2){
                    phosphate = true;
                }
            }

            if(carboxylicAcid){
                setType(24);
            }
            else if(phosphate){
                setType(24);
            }
            else if(imineOrEnol){
                setType(29);
            }
            else{
                setType(21);
            }
        }
        else if(neighborType == 7){
            setType(24);
        }
        else if(neighborType == 35){
            setType(21);
        }
        else if(neighborType == 49){
            setType(50);
        }
        else if(neighborType == 51){
            setType(52);
        }
        else if(neighborType == 70){
            setType(31);
        }
    }

    // phosphorus
    else if(neighbor->is(chemkit::Atom::Phosphorus)){
        setType(71);
    }

    // sulfur
    else if(neighbor->is(chemkit::Atom::Sulfur)){
        setType(71);
    }

    // silicon
    else if(neighbor->is(chemkit::Atom::Silicon)){
        setType(5);
    }
}

void MmffAtom::setCarbonType()
{
    const chemkit::Atom *atom = this->atom();

    // four neighbors
    if(atom->neighborCount() == 4){
        if(atom->isInRing(3)){
            setType(22); // carbon in three membered ring
        }
        else if(atom->isInRing(4)){
            if(atom->isBondedTo(chemkit::Atom::Carbon, chemkit::Bond::Double))
                setType(30); // olefinic carbon in four membered ring
            else
                setType(20); // carbon in foud membered ring
        }
        else{
            setType(1);
        }
    }

    // three neighbors
    else if(atom->neighborCount() == 3){
        const chemkit::Ring *smallestRing = atom->smallestRing();

        if(atom->isBondedTo(chemkit::Atom::Oxygen, chemkit::Bond::Double)){
            if(atom->neighborCount(chemkit::Atom::Oxygen) == 2){
                bool isNegative = false;
                foreach(const chemkit::Atom *neighbor, atom->neighbors()){
                    if(neighbor->is(chemkit::Atom::Oxygen) && neighbor->formalCharge() < 0){
                        isNegative = true;
                    }
                }

                if(isNegative){
                    setType(41);
                }
                else{
                    setType(3);
                }
            }
            else if(atom->neighborCount(chemkit::Atom::Nitrogen) == 1){
                setType(3); // amide carbonyl carbon
            }
            else if(atom->neighborCount(chemkit::Atom::Nitrogen) == 2){
                setType(3); // urea carbonyl carbon
            }
            else if(atom->neighborCount(chemkit::Atom::Carbon) >= 1){
                setType(3);
            }
            else{
                setType(3); // general carbonyl carbon
            }
        }
        else if(atom->isBondedTo(chemkit::Atom::Carbon, chemkit::Bond::Double)){
            if(atom->isInRing(4)){
                setType(30);
            }
            else{
                setType(2); // vinylic carbon
            }
        }
        else if(atom->isInRing() && smallestRing->size() == 3 && !smallestRing->isHeterocycle()){
            setType(22);
        }
        else if(isResonant(atom)){
            setType(57); // +N=C-N resonance structure
        }
        else if(isGuanidinium(atom)){
            setType(57); // CGD+ guanidinium
        }
        else if(atom->isBondedTo(chemkit::Atom::Nitrogen, chemkit::Bond::Double)){
            setType(3);
        }
        else if(smallestRing && smallestRing->size() == 4){
            setType(20);
        }
        else if(atom->isBondedTo(chemkit::Atom::Phosphorus, chemkit::Bond::Double) ||
                atom->isBondedTo(chemkit::Atom::Sulfur, chemkit::Bond::Double)){

            bool negativeSulfur = false;

            foreach(const chemkit::Atom *neighbor, atom->neighbors()){
                if(neighbor->is(chemkit::Atom::Sulfur) && neighbor->formalCharge() < 0)
                    negativeSulfur = true;
            }

            if(negativeSulfur && atom->neighborCount(chemkit::Atom::Sulfur) == 2)
                setType(41);
            else
                setType(3);
        }
        else{
            setType(2); // generic sp2 carbon
        }
    }

    // two neighbors
    else if(atom->neighborCount() == 2){
        if(atom->isBondedTo(chemkit::Atom::Nitrogen, 3) && atom->formalCharge() == -1){
            setType(60); // isonitrile carbon
        }
        else{
            setType(4); // acetylenic carbon
        }
    }

    // one neighbor
    else if(atom->neighborCount() == 1){
        if(atom->isBondedTo(chemkit::Atom::Nitrogen, 3) && atom->formalCharge() == -1){
            setType(60); // isonitrile carbon
        }
    }
}

void MmffAtom::setNitrogenType()
{
    const chemkit::Atom *atom = this->atom();

    // one neighbor
    if(atom->neighborCount() == 1){
        const chemkit::Bond *neighborBond = atom->bonds()[0];
        const chemkit::Atom *neighbor = neighborBond->otherAtom(atom);

        if(neighbor->is(chemkit::Atom::Carbon)){
            if(neighbor->isBondedTo(chemkit::Atom::Carbon, chemkit::Bond::Triple)){
                setType(40);
            }
            else if(neighbor->isBondedTo(chemkit::Atom::Nitrogen, chemkit::Bond::Double)){
                setType(40);
            }
            else{
                setType(42);
            }
        }
        else if(neighbor->is(chemkit::Atom::Nitrogen) && neighborBond->order() == chemkit::Bond::Double){
            setType(47);
        }
        else{
            setType(42);
        }
    }

    // two neighbors
    else if(atom->neighborCount() == 2){
        bool negativeRingNitrogen = false;

        if(atom->smallestRing() && atom->smallestRing()->size() == 5){
            foreach(const chemkit::Atom *ringAtom, atom->smallestRing()->atoms()){
                if(ringAtom->is(chemkit::Atom::Nitrogen) && ringAtom->formalCharge() == -1){
                    negativeRingNitrogen = true;
                }
            }
        }

        if(atom->isBondedTo(chemkit::Atom::Carbon, chemkit::Bond::Double) &&
           atom->isBondedTo(chemkit::Atom::Nitrogen, chemkit::Bond::Double)){
            setType(53);
        }
        else if(atom->formalCharge() == -1 || negativeRingNitrogen){
            setType(62, -1.0); // NM
        }
        else if(atom->isBondedTo(chemkit::Atom::Carbon, chemkit::Bond::Double)){
            setType(9);
        }
        else if(atom->isBondedTo(chemkit::Atom::Nitrogen, chemkit::Bond::Double)){
            int doubleBondedNitrogen = 0;

            foreach(const chemkit::Bond *bond, atom->bonds()){
                if(bond->otherAtom(atom)->is(chemkit::Atom::Nitrogen) &&
                   bond->order() == chemkit::Bond::Double){
                    doubleBondedNitrogen++;
                }
            }

            if(doubleBondedNitrogen == 2){
                setType(53);
            }
            else{
                setType(9);
            }
        }
        else if(atom->isInRing(5)){
            setType(79);
        }
        else if(atom->isBondedTo(chemkit::Atom::Oxygen, chemkit::Bond::Double)){
            setType(46); // nitroso
        }
        else if(atom->isBondedTo(chemkit::Atom::Carbon, chemkit::Bond::Triple)){
            setType(61); // isonitrile
        }
        else if(atom->isBondedTo(chemkit::Atom::Nitrogen, chemkit::Bond::Triple)){
            setType(61, 1.0); // diazo
        }
        else if(atom->isBondedTo(chemkit::Atom::Sulfur)){
            bool sulfate = false;
            bool nso = false;

            foreach(const chemkit::Atom *neighbor, atom->neighbors()){
                if(neighbor->is(chemkit::Atom::Sulfur) && isSulfate(neighbor)){
                    sulfate = true;
                }
                else if(neighbor->is(chemkit::Atom::Sulfur) &&
                        neighbor->isBondedTo(chemkit::Atom::Oxygen, chemkit::Bond::Double) &&
                        atom->bondTo(neighbor)->order() == chemkit::Bond::Double){
                    nso = true;
                }
            }

            if(sulfate){
                setType(43); // NSO2, NSO3
            }
            else if(nso){
                setType(48); // NSO
            }
            else{
                setType(8); // NR
            }
        }
        else{
            setType(8); // NR
        }
    }

    // three neighbors
    else if(atom->neighborCount() == 3){
        bool sulfate = false;
        bool phosphate = false;
        bool oxide = false;

        foreach(const chemkit::Atom *neighbor, atom->neighbors()){
            if(neighbor->is(chemkit::Atom::Sulfur) && isSulfate(neighbor)){
                sulfate = true;
            }
            else if(neighbor->is(chemkit::Atom::Phosphorus) && isPhosphate(neighbor)){
                phosphate = true;
            }
            else if(neighbor->is(chemkit::Atom::Oxygen) && neighbor->formalCharge() < 0){
                oxide = true;
            }
        }

        if(atom->isBondedTo(chemkit::Atom::Oxygen, chemkit::Bond::Double) && atom->neighborCount(chemkit::Atom::Oxygen) > 1){
            setType(45); // nitro or nitrate group nitrogen
        }
        else if(atom->formalCharge() == 1 && oxide){
            setType(67); // sp2 n-oxide nitrogen
        }
        else if(isGuanidinium(atom)){
            setType(56, (1.0/3.0)); // NGD+
        }
        else if(isResonant(atom)){
            setType(55, (1.0/2.0)); // NCN+
        }
        else if(atom->formalCharge() == 1 &&
                (atom->isBondedTo(chemkit::Atom::Carbon, chemkit::Bond::Double) ||
                 atom->isBondedTo(chemkit::Atom::Nitrogen, chemkit::Bond::Double))){
                setType(54, 1.0); // N+=C, N+=N
        }
        else if(sulfate || phosphate){
            setType(43);
        }
        else if(isAmide(atom)){
            setType(10);
        }
        else if(atom->isBondedTo(chemkit::Atom::Carbon)){
            bool doubleBond = false;
            bool doubleNitrogenBond = false;
            bool doubleNitrogenCarbonBond = false;
            bool cyano = false;

            foreach(const chemkit::Atom *neighbor, atom->neighbors()){
                if(neighbor->is(chemkit::Atom::Carbon)){
                    if(neighbor->isBondedTo(chemkit::Atom::Carbon, chemkit::Bond::Double) ||
                            neighbor->isBondedTo(chemkit::Atom::Nitrogen, chemkit::Bond::Double) ||
                            neighbor->isBondedTo(chemkit::Atom::Phosphorus, chemkit::Bond::Double)){
                        doubleBond = true;
                    }
                    else if(neighbor->isBondedTo(chemkit::Atom::Nitrogen, chemkit::Bond::Triple)){
                        cyano = true;
                    }
                }
                else if(neighbor->is(chemkit::Atom::Nitrogen)){
                    if(neighbor->isBondedTo(chemkit::Atom::Nitrogen, chemkit::Bond::Double)){
                        doubleNitrogenBond = true;
                    }
                    else if(neighbor->isBondedTo(chemkit::Atom::Carbon, chemkit::Bond::Double)){
                        doubleNitrogenCarbonBond = true;
                    }
                }
            }

            if(doubleBond){
                setType(40);
            }
            else if(doubleNitrogenBond){
                setType(10); // NN=N
            }
            else if(doubleNitrogenCarbonBond && !atom->isBondedTo(chemkit::Atom::Carbon)){
                setType(10); // NN=C
            }
            else if(cyano){
                setType(43); // nitrogen attached to cyano group
            }
            else{
                setType(8);
            }
        }
        else{
            setType(8); // nitrogen in aliphatic amines
        }
    }

    // four neighbors
    else if(atom->neighborCount() == 4){
        if(atom->neighborCount(chemkit::Atom::Oxygen) == 3){
            setType(45);
        }
        else if(atom->isBondedTo(chemkit::Atom::Oxygen, chemkit::Bond::Single)){
            foreach(const chemkit::Atom *neighbor, atom->neighbors()){
                if(neighbor->is(chemkit::Atom::Oxygen) && neighbor->formalCharge() == -1){
                    setType(68);
                }
            }
        }
        else{
            setType(34, 1.0);
        }
    }
}

void MmffAtom::setOxygenType()
{
    const chemkit::Atom *atom = this->atom();

    // one neighbor
    if(atom->neighborCount() == 1){
        const chemkit::Bond *neighborBond = atom->bonds()[0];
        const chemkit::Atom *neighbor = neighborBond->otherAtom(atom);

        if(neighbor->is(chemkit::Atom::Carbon)){
            if(neighborBond->order() == chemkit::Bond::Single){
                if(neighbor->isBondedTo(chemkit::Atom::Oxygen, chemkit::Bond::Double)){
                    if(atom->formalCharge() < 0){
                        setType(32, -0.5);
                    }
                    else{
                        setType(6);
                    }
                }
                else if(atom->formalCharge() < 0){
                    setType(35, -1.0); // alkoxide oxygen (OM)
                }
                else if(neighbor->isBondedTo(chemkit::Atom::Carbon, chemkit::Bond::Double)){
                    setType(6);
                }
                else if(neighbor->isBondedTo(chemkit::Atom::Nitrogen, chemkit::Bond::Double)){
                    setType(6);
                }
                else if(neighbor->isBondedTo(chemkit::Atom::Sulfur, chemkit::Bond::Double)){
                    setType(6);
                }
            }
            else if(neighborBond->order() == chemkit::Bond::Double){
                if(neighbor->isBondedTo(chemkit::Atom::Nitrogen)){
                    setType(7);
                }
                else if(neighbor->neighborCount(chemkit::Atom::Oxygen) > 1){
                    bool isNegative = false;
                    foreach(const chemkit::Atom *secondNeighbor, neighbor->neighbors()){
                        if(secondNeighbor->is(chemkit::Atom::Oxygen)){
                            if(secondNeighbor->formalCharge() < 0){
                                isNegative = true;
                            }
                        }
                    }
                    if(isNegative){
                        setType(32, -0.5);
                    }
                    else{
                        setType(7);
                    }
                }
                else if(neighbor->isBondedTo(chemkit::Atom::Carbon) || neighbor->isBondedTo(chemkit::Atom::Hydrogen)){
                    setType(7);
                }
                else{
                    setType(7);
                }
            }
        }
        else if(neighbor->is(chemkit::Atom::Nitrogen)){
            int oxygenCount = neighbor->neighborCount(chemkit::Atom::Oxygen);
            int negativeOxygenCount = 0;

            foreach(const chemkit::Atom *secondNeighbor, neighbor->neighbors()){
                if(secondNeighbor->is(chemkit::Atom::Oxygen) && secondNeighbor->formalCharge() < 0){
                    negativeOxygenCount++;
                }
            }

            if(oxygenCount >= 2){
                if(neighbor->is(chemkit::Atom::Nitrogen) && negativeOxygenCount == 1){
                    setType(32);
                }
                else if(neighbor->is(chemkit::Atom::Nitrogen) && oxygenCount == 3 && negativeOxygenCount > 1){
                    setType(32, -1.0/3.0);
                }
                else if(negativeOxygenCount > 1){
                    setType(32, -1.0/negativeOxygenCount);
                }
                else{
                    setType(32);
                }
            }
            else if(atom->formalCharge() < 0){
                if(neighbor->isBondedTo(chemkit::Atom::Carbon) &&
                   neighbor->neighborCount() == 2){
                    setType(35, -1.0);
                }
                else if(neighbor->formalCharge() == 0 && neighbor->neighborCount(chemkit::Atom::Oxygen) == 1){
                    setType(35, -1.0);
                }
                else{
                    setType(32);
                }
            }
            else{
                setType(7);
            }
        }
        else if(neighbor->is(chemkit::Atom::Sulfur)){
            int singleBondedOxygenCount = 0;
            int doubleBondedOxygenCount = 0;
            bool negativeOxygen = false;

            foreach(const chemkit::Atom *secondNeighbor, neighbor->neighbors()){
                if(secondNeighbor->is(chemkit::Atom::Oxygen)){
                    const chemkit::Bond *bond = neighbor->bondTo(secondNeighbor);

                    if(bond->order() == chemkit::Bond::Single){
                        singleBondedOxygenCount++;
                    }
                    else if(bond->order() == chemkit::Bond::Double){
                        doubleBondedOxygenCount++;
                    }

                    if(secondNeighbor->formalCharge() == -1){
                        negativeOxygen = true;
                    }
                }
            }

            int oxygenCount = singleBondedOxygenCount + doubleBondedOxygenCount;

            if(oxygenCount == 1 && neighbor->neighborCount() == 4){
                setType(32); // O-S
            }
            else if(doubleBondedOxygenCount >= 2){
                if(negativeOxygen){
                    setType(32, -1.0/3.0);
                }
                else if(neighbor->valence() == 5 && doubleBondedOxygenCount == 2){
                    setType(32, -0.5);
                }
                else{
                    setType(32); // O2S, O3S, 04S
                }
            }
            else if(singleBondedOxygenCount == 1 && doubleBondedOxygenCount == 1){
                setType(7);
            }
            else if(doubleBondedOxygenCount == 1 &&
                    neighbor->isBondedTo(chemkit::Atom::Sulfur, chemkit::Bond::Double) &&
                    neighbor->valence() == 5){
                setType(32, -0.5); // OSMS
            }
            else{
                setType(7);
            }
        }
        else if(neighbor->is(chemkit::Atom::Phosphorus)){
            int negativeOxygenAndSulfurCount = 0;
            bool doubleBondedOxygenOrSulfur = false;
            int oxygenAndSulfurCount = 0;

            foreach(const chemkit::Atom *secondNeighbor, neighbor->neighbors()){
                if(secondNeighbor->is(chemkit::Atom::Oxygen) || secondNeighbor->is(chemkit::Atom::Sulfur)){
                    oxygenAndSulfurCount++;

                    if(neighbor->bondTo(secondNeighbor)->order() == chemkit::Bond::Double){
                        doubleBondedOxygenOrSulfur = true;
                    }

                    if(secondNeighbor->isTerminal() && secondNeighbor->formalCharge() == -1){
                        negativeOxygenAndSulfurCount++;
                    }
                }
            }

            if(oxygenAndSulfurCount > 1 && doubleBondedOxygenOrSulfur && negativeOxygenAndSulfurCount){
                if(neighbor->valence() == 5 && negativeOxygenAndSulfurCount == 2){
                    setType(32, -2.0/3.0);
                }
                else{
                    setType(32, -0.5);
                }
            }
            else if(negativeOxygenAndSulfurCount > 1){
                setType(32, -1.0 / negativeOxygenAndSulfurCount);
            }
            else{
                setType(32);
            }
        }
        else if(neighbor->is(chemkit::Atom::Chlorine) && neighbor->neighborCount(chemkit::Atom::Oxygen) == 4){
            setType(32, -0.25); // O4CL
        }
        else if(neighbor->is(chemkit::Atom::Hydrogen) && atom->formalCharge() == -1){
            setType(35, -1.0);
        }
    }

    // two neighbors
    else if(atom->neighborCount() == 2){
        if(atom->neighborCount(chemkit::Atom::Hydrogen) == 2){
            setType(70); // water
        }
        else if(atom->formalCharge() == 1){
            if(atom->isBondedTo(chemkit::Atom::Carbon, chemkit::Bond::Double)){
                setType(51, 1.0);
            }
            else{
                setType(49, 1.0);
            }
        }
        else if(atom->isBondedTo(chemkit::Atom::Nitrogen)){
            setType(6);
        }
        else if(atom->isBondedTo(chemkit::Atom::Sulfur)){
            setType(6);
        }
        else if(atom->isBondedTo(chemkit::Atom::Phosphorus)){
            setType(6);
        }
        else if(atom->isBondedTo(chemkit::Atom::Carbon)){
            setType(6);
        }
        else{
            setType(6);
        }
    }

    // three neighbors
    else if(atom->neighborCount() == 3){
        setType(49, 1.0);
    }
}

void MmffAtom::setSulfurType()
{
    const chemkit::Atom *atom = this->atom();

    if(atom->isTerminal()){
        const chemkit::Atom *neighbor = atom->bonds()[0]->otherAtom(atom);

        if(isThiocarboxylate(atom)){
            setType(72, -0.5);
        }
        else if(atom->isBondedTo(chemkit::Atom::Carbon, chemkit::Bond::Double)){
            setType(16);
        }
        else if(atom->isBondedTo(chemkit::Atom::Oxygen, chemkit::Bond::Double)){
            setType(17);
        }
        else if(neighbor->is(chemkit::Atom::Phosphorus)){
            if(neighbor->isBondedTo(chemkit::Atom::Oxygen, chemkit::Bond::Double) && atom->formalCharge() == -1){
                setType(72, -0.5);
            }
            else{
                setType(72); // S-P
            }
        }
        else if(atom->formalCharge() < 0){
            setType(72, -1.0); // SM
        }
        else if(atom->isBondedTo(chemkit::Atom::Sulfur, chemkit::Bond::Double)){
            setType(72, -0.5);
        }
        else{
            setType(72);
        }
    }
    else if(atom->isBondedTo(chemkit::Atom::Nitrogen, chemkit::Bond::Double) &&
            atom->isBondedTo(chemkit::Atom::Oxygen, chemkit::Bond::Double)){
        setType(18);
    }
    else if(atom->isBondedTo(chemkit::Atom::Nitrogen, chemkit::Bond::Double) && atom->neighborCount() == 3){
        setType(17); // >S=N
    }
    else if(atom->isBondedTo(chemkit::Atom::Oxygen, chemkit::Bond::Double) &&
            atom->isBondedTo(chemkit::Atom::Sulfur, chemkit::Bond::Double)){
        setType(73);
    }
    else if(atom->isBondedTo(chemkit::Atom::Oxygen, chemkit::Bond::Double)){
        int singleBondedOxygenCount = 0;
        int doubleBondedOxygenCount = 0;

        foreach(const chemkit::Atom *neighbor, atom->neighbors()){
            if(neighbor->is(chemkit::Atom::Oxygen)){
                const chemkit::Bond *bond = atom->bondTo(neighbor);

                if(bond->order() == chemkit::Bond::Single){
                    singleBondedOxygenCount++;
                }
                else if(bond->order() == chemkit::Bond::Double){
                    doubleBondedOxygenCount++;
                }
            }
        }

        if(singleBondedOxygenCount == 1 && doubleBondedOxygenCount == 1){
            setType(17); // S=O
        }
        else if(doubleBondedOxygenCount == 2 && atom->valence() == 5){
            setType(73); // SO2M
        }
        else if(doubleBondedOxygenCount == 1 && atom->isBondedTo(chemkit::Atom::Carbon, chemkit::Bond::Double)){
            setType(74); // =S=O
        }
        else if(doubleBondedOxygenCount >= 2){
            setType(18); // SO2, SO2N, SO3, SO4, =SO2, SNO
        }
        else{
            setType(17); // S=O
        }
    }
    else{
        setType(15);
    }
}

void MmffAtom::setAromaticType(const chemkit::Ring *aromaticRing, int position)
{
    const chemkit::Atom *atom = this->atom();
    int type = this->typeNumber();

    // carbon
    if(atom->is(chemkit::Atom::Carbon)){
        if(aromaticRing->size() == 5){
            if(type == 57){
                setType(80); // CIM+
            }
            else if(position == 0){
                setType(78); // C5
            }
            else if(position == 1){
                if(type == 64){
                    setType(78); // C5
                }
                else{
                    setType(63); // C5A
                }
            }
            else if(position == 2){
                if(type == 63){
                    setType(78); // C5
                }
                else{
                    setType(64); // C5B
                }
            }
            else{
                setType(78); // C5
            }
        }
        else if(aromaticRing->size() == 6){
            setType(37); // CB
        }
    }

    // nitrogen
    else if(atom->is(chemkit::Atom::Nitrogen)){
        if(aromaticRing->size() == 5){
            if(type == 62){
                if(aromaticRing->atomCount(chemkit::Atom::Nitrogen) == 2){
                    setType(76, -0.5); // N5M
                }
                else if(aromaticRing->atomCount(chemkit::Atom::Nitrogen) == 3){
                    setType(76, -1.0/3.0); // N5M
                }
                else if(aromaticRing->atomCount(chemkit::Atom::Nitrogen) == 4){
                    setType(76, -1.0/4.0); // N5M
                }
            }
            else if(type == 67){
                setType(82); // N5OX
            }
            else if(type == 54){
                setType(81, 1.0); // N5+
            }
            else if(type == 55){
                setType(81, 0.5); // NIM+
            }
            else if(type == 56){
                setType(81, 1.0/3.0);
            }
            else if(position == 0){
                setType(39); // NPYL
            }
            else if(position == 1){
                if(type == 66){
                    setType(79); // N5
                }
                else{
                    setType(65); // N5A
                }
            }
            else if(position == 2){
                if(type == 65){
                    setType(79); // N5
                }
                else{
                    setType(66); // N5A
                }
            }
            else{
                setType(79); // N5
            }
        }
        else if(aromaticRing->size() == 6){
            if(type == 54 || type == 55 || type == 56){
                setType(58, 1.0); // NPD+
            }
            else if(type == 67){
                setType(69); // NPOX
            }
            else{
                if(atom->formalCharge() > 0){
                    setType(58, 1.0); // NPYD+
                }
                else{
                    setType(38); // NPYD
                }
            }
        }
    }

    // oxygen
    else if(atom->is(chemkit::Atom::Oxygen)){
        if(aromaticRing->size() == 5){
            setType(59); // OFUR
        }
    }

    // sulfur
    else if(atom->is(chemkit::Atom::Sulfur)){
        if(aromaticRing->size() == 5){
            setType(44); // STHI
        }
    }
}

// --- Parameters ---------------------------------------------------------- //
const MmffAtomParameters* MmffAtom::parameters() const
{
    const MmffParameters *mmffParameters = static_cast<const MmffForceField *>(forceField())->parameters();
    if(!mmffParameters)
        return 0;

    return mmffParameters->atomParameters(this);
}

// --- Internal Methods ---------------------------------------------------- //
bool MmffAtom::isGuanidinium(const chemkit::Atom *atom) const
{
    if(atom->is(chemkit::Atom::Carbon)){
        bool doubleBondedPositiveNitrogen = false;
        int singleBondedNitrogenCount = 0;

        foreach(const chemkit::Atom *neighbor, atom->neighbors()){
            if(neighbor->is(chemkit::Atom::Nitrogen)){
                if(neighbor->formalCharge() == 1 &&
                   atom->bondTo(neighbor)->order() == chemkit::Bond::Double){
                    doubleBondedPositiveNitrogen = true;
                }
                else if(atom->bondTo(neighbor)->order() == chemkit::Bond::Single &&
                        neighbor->neighborCount() == 3){
                    singleBondedNitrogenCount++;
                }
            }
        }

        if(doubleBondedPositiveNitrogen && singleBondedNitrogenCount == 2){
            return true;
        }
    }
    else if(atom->is(chemkit::Atom::Nitrogen)){
        foreach(const chemkit::Atom *neighbor, atom->neighbors()){
            if(neighbor->is(chemkit::Atom::Carbon)){
                if(isGuanidinium(neighbor)){
                    return true;
                }
            }
        }
    }

    return false;
}

bool MmffAtom::isResonant(const chemkit::Atom *atom) const
{
    if(atom->is(chemkit::Atom::Carbon)){
        const chemkit::Atom *doubleBondedPositiveNitrogen = 0;
        int singleBondedNitrogenCount = 0;

        foreach(const chemkit::Atom *neighbor, atom->neighbors()){
            const chemkit::Bond *bond = atom->bondTo(neighbor);

            if(neighbor->is(chemkit::Atom::Nitrogen)){
                if(bond->order() == chemkit::Bond::Double &&
                   neighbor->neighborCount() == 3 &&
                   neighbor->formalCharge() == 1){
                    doubleBondedPositiveNitrogen = neighbor;

                    foreach(const chemkit::Atom *secondNeighbor, neighbor->neighbors()){
                        if(secondNeighbor->formalCharge() < 0 && secondNeighbor->isTerminal()){
                            doubleBondedPositiveNitrogen = 0;
                            break;
                        }
                    }
                }
                else if(bond->order() == chemkit::Bond::Single &&
                        neighbor->neighborCount() == 3 &&
                        neighbor->formalCharge() == 0){
                    singleBondedNitrogenCount++;
                }
            }
        }

        foreach(const chemkit::Ring *ring, atom->rings()){
            if(ring->contains(doubleBondedPositiveNitrogen) && isPositiveAromaticNitrogenRing(ring)){
                return false;
            }
        }

        if(doubleBondedPositiveNitrogen && singleBondedNitrogenCount == 1){
            return true;
        }
    }
    else if(atom->is(chemkit::Atom::Nitrogen)){
        foreach(const chemkit::Atom *neighbor, atom->neighbors()){
            if(neighbor->is(chemkit::Atom::Carbon)){
                if(isResonant(neighbor)){
                    return true;
                }
            }
        }
    }

    return false;
}

bool MmffAtom::isAmide(const chemkit::Atom *atom) const
{
    if(atom->is(chemkit::Atom::Carbon)){
        if(atom->isBondedTo(chemkit::Atom::Oxygen, chemkit::Bond::Double) &&
           atom->isBondedTo(chemkit::Atom::Nitrogen, chemkit::Bond::Single)){
            return true;
        }
        else if(atom->isBondedTo(chemkit::Atom::Sulfur, chemkit::Bond::Double) &&
                atom->isBondedTo(chemkit::Atom::Nitrogen, chemkit::Bond::Single)){
            return true;
        }
    }
    else if(atom->is(chemkit::Atom::Nitrogen)){
        foreach(const chemkit::Atom *neighbor, atom->neighbors()){
            if(neighbor->is(chemkit::Atom::Carbon)){
                if(isAmide(neighbor)){
                    return true;
                }
            }
        }
    }

    return false;
}

bool MmffAtom::isPhosphate(const chemkit::Atom *atom) const
{
    if(atom->is(chemkit::Atom::Phosphorus)){
        if(atom->isInRing())
            return false;

        int singleBondedOxygenCount = 0;
        int doubleBondedOxygenCount = 0;

        foreach(const chemkit::Atom *neighbor, atom->neighbors()){
            if(neighbor->is(chemkit::Atom::Oxygen)){
                const chemkit::Bond *bond = atom->bondTo(neighbor);

                if(bond->order() == chemkit::Bond::Single){
                    singleBondedOxygenCount++;
                }
                else if(bond->order() == chemkit::Bond::Double){
                    doubleBondedOxygenCount++;
                }
            }
        }

        int oxygenCount = singleBondedOxygenCount + doubleBondedOxygenCount;
        if(oxygenCount >= 2 &&
           doubleBondedOxygenCount >= 1){
            return true;
        }
    }

    return false;
}

bool MmffAtom::isSulfate(const chemkit::Atom *atom) const
{
    if(atom->is(chemkit::Atom::Sulfur)){
        int singleBondedOxygenCount = 0;
        int doubleBondedOxygenCount = 0;

        foreach(const chemkit::Atom *neighbor, atom->neighbors()){
            if(neighbor->is(chemkit::Atom::Oxygen)){
                const chemkit::Bond *bond = atom->bondTo(neighbor);

                if(bond->order() == chemkit::Bond::Single){
                    singleBondedOxygenCount++;
                }
                else if(bond->order() == chemkit::Bond::Double){
                    doubleBondedOxygenCount++;
                }
            }
        }

        int oxygenCount = singleBondedOxygenCount + doubleBondedOxygenCount;
        if(oxygenCount >= 2 && doubleBondedOxygenCount >= 1){
            return true;
        }
    }

    return false;
}

bool MmffAtom::isThiocarboxylate(const chemkit::Atom *atom) const
{
    if(atom->is(chemkit::Atom::Carbon)){
        bool negativeSulfur = false;
        bool doubleBondedSulfur = false;
        int sulfurCount = 0;

        foreach(const chemkit::Atom *neighbor, atom->neighbors()){
            const chemkit::Bond *bond = atom->bondTo(neighbor);

            if(neighbor->is(chemkit::Atom::Sulfur) && neighbor->isTerminal()){
                sulfurCount++;

                if(bond->order() == chemkit::Bond::Single && neighbor->formalCharge() == -1){
                    negativeSulfur = true;
                }
                else if(bond->order() == chemkit::Bond::Double && neighbor->formalCharge() == 0){
                    doubleBondedSulfur = true;
                }
            }
        }

        if(sulfurCount == 2 && negativeSulfur && doubleBondedSulfur){
            return true;
        }
    }
    else if(atom->is(chemkit::Atom::Sulfur)){
        foreach(const chemkit::Atom *neighbor, atom->neighbors()){
            if(neighbor->is(chemkit::Atom::Carbon) && isThiocarboxylate(neighbor)){
                return true;
            }
        }
    }

    return false;
}

bool MmffAtom::isPositiveAromaticNitrogenRing(const chemkit::Ring *ring) const
{
    if(ring->size() != 6){
        return false;
    }

    int nitrogenCount = 0;
    const chemkit::Atom *positiveNitrogen = 0;
    const chemkit::Atom *neutralNitrogen = 0;

    foreach(const chemkit::Atom *atom, ring->atoms()){
        if(atom->is(chemkit::Atom::Nitrogen)){
            nitrogenCount++;

            if(atom->formalCharge() == 1 &&
               atom->neighborCount() == 3){
                positiveNitrogen = atom;
            }
            else if(atom->formalCharge() == 0 &&
                    atom->neighborCount() == 2){
                neutralNitrogen = atom;
            }
        }
    }

    if(!positiveNitrogen){
        return false;
    }

    int doubleBondCount = 0;
    foreach(const chemkit::Bond *bond, ring->bonds()){
        if(bond->order() == chemkit::Bond::Double){
            doubleBondCount++;
        }
    }

    if(doubleBondCount != 3){
        return false;
    }

    return true;
}
