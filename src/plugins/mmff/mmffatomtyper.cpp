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

#include "mmffatomtyper.h"

#include <boost/math/special_functions/round.hpp>

#include <chemkit/atom.h>
#include <chemkit/foreach.h>
#include <chemkit/molecule.h>

#include "mmffaromaticitymodel.h"

namespace {

bool isGuanidinium(const chemkit::Atom *atom)
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

bool isPositiveAromaticNitrogenRing(const chemkit::Ring *ring)
{
    if(ring->size() != 6){
        return false;
    }

    int nitrogenCount = 0;
    const chemkit::Atom *positiveNitrogen = 0;

    foreach(const chemkit::Atom *atom, ring->atoms()){
        if(atom->is(chemkit::Atom::Nitrogen)){
            nitrogenCount++;

            if(atom->formalCharge() == 1 &&
               atom->neighborCount() == 3){
                positiveNitrogen = atom;
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

bool isResonant(const chemkit::Atom *atom)
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

bool isAmide(const chemkit::Atom *atom)
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

bool isPhosphate(const chemkit::Atom *atom)
{
    if(atom->is(chemkit::Atom::Phosphorus)){
        if(atom->isInRing()){
            return false;
        }

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

bool isSulfate(const chemkit::Atom *atom)
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

bool isThiocarboxylate(const chemkit::Atom *atom)
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

int ringPosition(const chemkit::Atom *atom, const chemkit::Ring *ring)
{
    if(ring->size() != 5){
        return 0;
    }

    const chemkit::Atom *ringRoot = 0;
    int rootAtomCount = 0;

    bool imidizadole = false;
    bool positiveNitrogen = false;

    if(ring->size() == 5){
        foreach(const chemkit::Atom *atom, ring->atoms()){
            if(atom->is(chemkit::Atom::Nitrogen) &&
               atom->neighborCount() == 3 &&
               atom->valence() == 3){

                if(ringRoot && ringRoot->atomicNumber() == atom->atomicNumber()){
                    rootAtomCount++;
                }
                else if(ringRoot && atom->atomicNumber() > ringRoot->atomicNumber()){
                    ringRoot = atom;
                    rootAtomCount++;
                }
                else if(ringRoot && atom->atomicNumber() < ringRoot->atomicNumber()){
                }
                else{
                    ringRoot = atom;
                    rootAtomCount = 0;
                }
            }
            else if(atom->is(chemkit::Atom::Nitrogen) && atom->formalCharge() == 1){
                bool negativeNeighbor = false;

                foreach(const chemkit::Atom *neighbor, atom->neighbors()){
                    if(neighbor->formalCharge() < 0){
                        negativeNeighbor = true;
                    }
                }

                if(!negativeNeighbor){
                    positiveNitrogen = true;
                }
            }
            else if((atom->is(chemkit::Atom::Oxygen) || atom->is(chemkit::Atom::Sulfur)) &&
                    atom->neighborCount() == 2){

                if(ringRoot && ringRoot->atomicNumber() == atom->atomicNumber()){
                    rootAtomCount++;
                }
                else if(ringRoot && atom->atomicNumber() > ringRoot->atomicNumber()){
                    ringRoot = atom;
                    rootAtomCount++;
                }
                else if(ringRoot && atom->atomicNumber() < ringRoot->atomicNumber()){
                }
                else{
                    ringRoot = atom;
                    rootAtomCount = 0;
                }
            }
        }
    }

    if(positiveNitrogen && ring->atomCount(chemkit::Atom::Nitrogen) >= 2){
        imidizadole = true;

        foreach(const chemkit::Atom *atom, ring->atoms()){
            if(atom->is(chemkit::Atom::Nitrogen) && atom->formalCharge() == 1){
                if(atom->isBondedTo(chemkit::Atom::Nitrogen)){
                    imidizadole = false;
                }
            }
        }
    }

    if(imidizadole && ring->heteroatomCount() == 2){
        return 0;
    }
    else if(rootAtomCount > 1){
        ringRoot = 0;
    }
    else if(!ringRoot){
        return 0;
    }
    else if(positiveNitrogen && ringRoot->is(chemkit::Atom::Nitrogen)){
        return 4;
    }

    return ring->position(atom, ringRoot);
}

} // end anonymous namespace

// --- Construction and Destruction ---------------------------------------- //
MmffAtomTyper::MmffAtomTyper(const chemkit::Molecule *molecule)
    : chemkit::AtomTyper("mmff")
{
    setMolecule(molecule);
}

MmffAtomTyper::~MmffAtomTyper()
{
}

// --- Properties ---------------------------------------------------------- //
void MmffAtomTyper::setMolecule(const chemkit::Molecule *molecule)
{
    chemkit::AtomTyper::setMolecule(molecule);

    if(!molecule){
        m_types.resize(0);
        return;
    }

    m_types.resize(molecule->atomCount());
    m_formalCharges.resize(molecule->atomCount());

    MmffAromaticityModel aromaticityModel;
    aromaticityModel.setMolecule(molecule);

    // assign types to heavy atoms
    foreach(const chemkit::Atom *atom, molecule->atoms()){
        if(atom->isTerminalHydrogen()){
            continue;
        }

        setType(atom->index(), atom);
    }

    // assign aromatic atom types
    std::vector<const chemkit::Ring *> sixMemberedAromaticRings;
    std::vector<const chemkit::Ring *> fiveMemberedAromaticRings;
    foreach(const chemkit::Ring *ring, molecule->rings()){
        if(ring->size() == 5 && aromaticityModel.isAromatic(ring)){
            fiveMemberedAromaticRings.push_back(ring);
        }
        else if(ring->size() == 6 && aromaticityModel.isAromatic(ring)){
            sixMemberedAromaticRings.push_back(ring);
        }
    }

    foreach(const chemkit::Ring *ring, sixMemberedAromaticRings){
        foreach(const chemkit::Atom *atom, ring->atoms()){
            setAromaticType(atom->index(), atom, ring, ringPosition(atom, ring));
        }
    }

    foreach(const chemkit::Ring *ring, fiveMemberedAromaticRings){
        foreach(const chemkit::Atom *atom, ring->atoms()){
            setAromaticType(atom->index(), atom, ring, ringPosition(atom, ring));
        }
    }

    // assign terminal hydrogen types
    foreach(const chemkit::Atom *atom, molecule->atoms()){
        if(atom->isTerminalHydrogen()){
            setHydrogenType(atom->index(), atom);
        }
    }
}

// --- Types --------------------------------------------------------------- //
int MmffAtomTyper::typeNumber(const chemkit::Atom *atom) const
{
    return m_types[atom->index()];
}

// --- Charges ------------------------------------------------------------- //
chemkit::Real MmffAtomTyper::formalCharge(int index) const
{
    return m_formalCharges[index];
}

chemkit::Real MmffAtomTyper::formalCharge(const chemkit::Atom *atom) const
{
    return formalCharge(atom->index());
}

// --- Internal Methods ---------------------------------------------------- //
void MmffAtomTyper::setType(int index, int type, chemkit::Real formalCharge)
{
    m_types[index] = type;
    m_formalCharges[index] = formalCharge;
}

void MmffAtomTyper::setType(int index, const chemkit::Atom *atom)
{
    switch(atom->atomicNumber()){
        // carbon
        case chemkit::Atom::Carbon:
            setCarbonType(index, atom);
            break;

        // nitrogen
        case chemkit::Atom::Nitrogen:
            setNitrogenType(index, atom);
            break;

        // oxygen
        case chemkit::Atom::Oxygen:
            setOxygenType(index, atom);
            break;

        // phosphorus
        case chemkit::Atom::Phosphorus:
            if(atom->neighborCount() == 4){
                setType(index, 25);
            }
            else if(atom->neighborCount() == 3){
                if(atom->isBondedTo(chemkit::Atom::Carbon, chemkit::Bond::Double)){
                    setType(index, 75);
                }
                else{
                    setType(index, 26);
                }
            }
            else if(atom->neighborCount() == 2 && atom->isBondedTo(chemkit::Atom::Carbon)){
                setType(index, 75);
            }
            break;

        // sulfur
        case chemkit::Atom::Sulfur:
            setSulfurType(index, atom);
            break;

        // fluorine
        case chemkit::Atom::Fluorine:
            if(atom->valence() > 0){
                setType(index, 11);
            }
            else{
                setType(index, 89, -1.0);
            }
            break;

        // chlorine
        case chemkit::Atom::Chlorine:
            if(atom->neighborCount(chemkit::Atom::Oxygen) == 4){
                setType(index, 77);
            }
            else if(atom->valence() > 0){
                setType(index, 12);
            }
            else{
                setType(index, 90, -1.0);
            }
            break;

        // bromine
        case chemkit::Atom::Bromine:
            if(atom->valence() > 0){
                setType(index, 13);
            }
            else{
                setType(index, 91, -1.0);
            }
            break;

        // iodine
        case chemkit::Atom::Iodine:
            if(atom->valence() > 0){
                setType(index, 14);
            }
            break;

        // iron
        case chemkit::Atom::Iron:
            if(boost::math::iround(atom->partialCharge()) == 2){
                setType(index, 87, 2.0);
            }
            else{
                setType(index, 88, 3.0);
            }
            break;

        // lithium
        case chemkit::Atom::Lithium:
            setType(index, 92, 1.0);
            break;

        // sodium
        case chemkit::Atom::Sodium:
            setType(index, 93, 1.0);
            break;

        // potassium
        case chemkit::Atom::Potassium:
            setType(index, 94, 1.0);
            break;

        // zinc
        case chemkit::Atom::Zinc:
            setType(index, 95, 2.0);
            break;

        // calcium
        case chemkit::Atom::Calcium:
            setType(index, 96, 2.0);
            break;

        // copper
        case chemkit::Atom::Copper:
            if(boost::math::iround(atom->partialCharge()) == 2){
                setType(index, 98, 2.0);
            }
            else{
                setType(index, 97, 1.0);
            }
            break;

        // magnesium
        case chemkit::Atom::Magnesium:
            setType(index, 99, 2.0);
            break;

        // silicon
        case chemkit::Atom::Silicon:
            setType(index, 19);
            break;

        default:
            break;
    }
}

void MmffAtomTyper::setHydrogenType(int index, const chemkit::Atom *atom)
{
    assert(atom->isTerminalHydrogen());

    const chemkit::Atom *neighbor = atom->neighbor(0);
    int neighborType = typeNumber(neighbor);

    // carbon
    if(neighbor->is(chemkit::Atom::Carbon)){
        setType(index, 5);
    }

    // nitrogen
    else if(neighbor->is(chemkit::Atom::Nitrogen)){
        if(neighborType == 8 || neighborType == 39 || neighborType == 45 ||
           neighborType == 62 || neighborType == 67 || neighborType == 68){
            setType(index, 23);
        }
        else if(neighborType == 9){
            setType(index, 27);
        }
        else if(neighborType == 10 || neighborType == 40 || neighborType == 42 ||
                neighborType == 43 || neighborType == 48){
            setType(index, 28);
        }
        else if(neighborType == 34 || neighborType == 54 || neighborType == 55 ||
                neighborType == 56 || neighborType == 58 || neighborType == 81){
            setType(index, 36);
        }
    }

    // oxygen
    else if(neighbor->is(chemkit::Atom::Oxygen)){
        if(neighbor->isBondedTo(chemkit::Atom::Sulfur)){
            setType(index, 33);
        }
        else if(neighborType == 6){
            bool imineOrEnol = false;
            bool carboxylicAcid = false;
            bool phosphate = false;

            foreach(const chemkit::Atom *secondNeighbor, neighbor->neighbors()){
                if(secondNeighbor == atom){
                    continue;
                }

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
                setType(index, 24);
            }
            else if(phosphate){
                setType(index, 24);
            }
            else if(imineOrEnol){
                setType(index, 29);
            }
            else{
                setType(index, 21);
            }
        }
        else if(neighborType == 7){
            setType(index, 24);
        }
        else if(neighborType == 35){
            setType(index, 21);
        }
        else if(neighborType == 49){
            setType(index, 50);
        }
        else if(neighborType == 51){
            setType(index, 52);
        }
        else if(neighborType == 70){
            setType(index, 31);
        }
    }

    // phosphorus
    else if(neighbor->is(chemkit::Atom::Phosphorus)){
        setType(index, 71);
    }

    // sulfur
    else if(neighbor->is(chemkit::Atom::Sulfur)){
        setType(index, 71);
    }

    // silicon
    else if(neighbor->is(chemkit::Atom::Silicon)){
        setType(index, 5);
    }
}

void MmffAtomTyper::setCarbonType(int index, const chemkit::Atom *atom)
{
    // four neighbors
    if(atom->neighborCount() == 4){
        if(atom->isInRing(3)){
            setType(index, 22); // carbon in three membered ring
        }
        else if(atom->isInRing(4)){
            if(atom->isBondedTo(chemkit::Atom::Carbon, chemkit::Bond::Double)){
                setType(index, 30); // olefinic carbon in four membered ring
            }
            else{
                setType(index, 20); // carbon in foud membered ring
            }
        }
        else{
            setType(index, 1);
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
                    setType(index, 41);
                }
                else{
                    setType(index, 3);
                }
            }
            else if(atom->neighborCount(chemkit::Atom::Nitrogen) == 1){
                setType(index, 3); // amide carbonyl carbon
            }
            else if(atom->neighborCount(chemkit::Atom::Nitrogen) == 2){
                setType(index, 3); // urea carbonyl carbon
            }
            else if(atom->neighborCount(chemkit::Atom::Carbon) >= 1){
                setType(index, 3);
            }
            else{
                setType(index, 3); // general carbonyl carbon
            }
        }
        else if(atom->isBondedTo(chemkit::Atom::Carbon, chemkit::Bond::Double)){
            if(atom->isInRing(4)){
                setType(index, 30);
            }
            else{
                setType(index, 2); // vinylic carbon
            }
        }
        else if(atom->isInRing() && smallestRing->size() == 3 && !smallestRing->isHeterocycle()){
            setType(index, 22);
        }
        else if(isResonant(atom)){
            setType(index, 57); // +N=C-N resonance structure
        }
        else if(isGuanidinium(atom)){
            setType(index, 57); // CGD+ guanidinium
        }
        else if(atom->isBondedTo(chemkit::Atom::Nitrogen, chemkit::Bond::Double)){
            setType(index, 3);
        }
        else if(smallestRing && smallestRing->size() == 4){
            setType(index, 20);
        }
        else if(atom->isBondedTo(chemkit::Atom::Phosphorus, chemkit::Bond::Double) ||
                atom->isBondedTo(chemkit::Atom::Sulfur, chemkit::Bond::Double)){

            bool negativeSulfur = false;

            foreach(const chemkit::Atom *neighbor, atom->neighbors()){
                if(neighbor->is(chemkit::Atom::Sulfur) && neighbor->formalCharge() < 0){
                    negativeSulfur = true;
                }
            }

            if(negativeSulfur && atom->neighborCount(chemkit::Atom::Sulfur) == 2){
                setType(index, 41);
            }
            else{
                setType(index, 3);
            }
        }
        else{
            setType(index, 2); // generic sp2 carbon
        }
    }

    // two neighbors
    else if(atom->neighborCount() == 2){
        if(atom->isBondedTo(chemkit::Atom::Nitrogen, 3) && atom->formalCharge() == -1){
            setType(index, 60); // isonitrile carbon
        }
        else{
            setType(index, 4); // acetylenic carbon
        }
    }

    // one neighbor
    else if(atom->neighborCount() == 1){
        if(atom->isBondedTo(chemkit::Atom::Nitrogen, 3) && atom->formalCharge() == -1){
            setType(index, 60); // isonitrile carbon
        }
    }
}

void MmffAtomTyper::setNitrogenType(int index, const chemkit::Atom *atom)
{
    // one neighbor
    if(atom->neighborCount() == 1){
        const chemkit::Bond *neighborBond = atom->bonds()[0];
        const chemkit::Atom *neighbor = neighborBond->otherAtom(atom);

        if(neighbor->is(chemkit::Atom::Carbon)){
            if(neighbor->isBondedTo(chemkit::Atom::Carbon, chemkit::Bond::Triple)){
                setType(index, 40);
            }
            else if(neighbor->isBondedTo(chemkit::Atom::Nitrogen, chemkit::Bond::Double)){
                setType(index, 40);
            }
            else{
                setType(index, 42);
            }
        }
        else if(neighbor->is(chemkit::Atom::Nitrogen) && neighborBond->order() == chemkit::Bond::Double){
            setType(index, 47);
        }
        else{
            setType(index, 42);
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
            setType(index, 53);
        }
        else if(atom->formalCharge() == -1 || negativeRingNitrogen){
            setType(index, 62, -1.0); // NM
        }
        else if(atom->isBondedTo(chemkit::Atom::Carbon, chemkit::Bond::Double)){
            setType(index, 9);
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
                setType(index, 53);
            }
            else{
                setType(index, 9);
            }
        }
        else if(atom->isInRing(5)){
            setType(index, 79);
        }
        else if(atom->isBondedTo(chemkit::Atom::Oxygen, chemkit::Bond::Double)){
            setType(index, 46); // nitroso
        }
        else if(atom->isBondedTo(chemkit::Atom::Carbon, chemkit::Bond::Triple)){
            setType(index, 61); // isonitrile
        }
        else if(atom->isBondedTo(chemkit::Atom::Nitrogen, chemkit::Bond::Triple)){
            setType(index, 61, 1.0); // diazo
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
                setType(index, 43); // NSO2, NSO3
            }
            else if(nso){
                setType(index, 48); // NSO
            }
            else{
                setType(index, 8); // NR
            }
        }
        else{
            setType(index, 8); // NR
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
            setType(index, 45); // nitro or nitrate group nitrogen
        }
        else if(atom->formalCharge() == 1 && oxide){
            setType(index, 67); // sp2 n-oxide nitrogen
        }
        else if(isGuanidinium(atom)){
            setType(index, 56, (1.0/3.0)); // NGD+
        }
        else if(isResonant(atom)){
            setType(index, 55, (1.0/2.0)); // NCN+
        }
        else if(atom->formalCharge() == 1 &&
                (atom->isBondedTo(chemkit::Atom::Carbon, chemkit::Bond::Double) ||
                 atom->isBondedTo(chemkit::Atom::Nitrogen, chemkit::Bond::Double))){
            setType(index, 54, 1.0); // N+=C, N+=N
        }
        else if(sulfate || phosphate){
            setType(index, 43);
        }
        else if(isAmide(atom)){
            setType(index, 10);
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
                setType(index, 40);
            }
            else if(doubleNitrogenBond){
                setType(index, 10); // NN=N
            }
            else if(doubleNitrogenCarbonBond && !atom->isBondedTo(chemkit::Atom::Carbon)){
                setType(index, 10); // NN=C
            }
            else if(cyano){
                setType(index, 43); // nitrogen attached to cyano group
            }
            else{
                setType(index, 8);
            }
        }
        else{
            setType(index, 8); // nitrogen in aliphatic amines
        }
    }

    // four neighbors
    else if(atom->neighborCount() == 4){
        if(atom->neighborCount(chemkit::Atom::Oxygen) == 3){
            setType(index, 45);
        }
        else if(atom->isBondedTo(chemkit::Atom::Oxygen, chemkit::Bond::Single)){
            foreach(const chemkit::Atom *neighbor, atom->neighbors()){
                if(neighbor->is(chemkit::Atom::Oxygen) && neighbor->formalCharge() == -1){
                    setType(index, 68);
                }
            }
        }
        else{
            setType(index, 34, 1.0);
        }
    }
}

void MmffAtomTyper::setOxygenType(int index, const chemkit::Atom *atom)
{
    // one neighbor
    if(atom->neighborCount() == 1){
        const chemkit::Bond *neighborBond = atom->bonds()[0];
        const chemkit::Atom *neighbor = neighborBond->otherAtom(atom);

        if(neighbor->is(chemkit::Atom::Carbon)){
            if(neighborBond->order() == chemkit::Bond::Single){
                if(neighbor->isBondedTo(chemkit::Atom::Oxygen, chemkit::Bond::Double)){
                    if(atom->formalCharge() < 0){
                        setType(index, 32, -0.5);
                    }
                    else{
                        setType(index, 6);
                    }
                }
                else if(atom->formalCharge() < 0){
                    setType(index, 35, -1.0); // alkoxide oxygen (OM)
                }
                else if(neighbor->isBondedTo(chemkit::Atom::Carbon, chemkit::Bond::Double)){
                    setType(index, 6);
                }
                else if(neighbor->isBondedTo(chemkit::Atom::Nitrogen, chemkit::Bond::Double)){
                    setType(index, 6);
                }
                else if(neighbor->isBondedTo(chemkit::Atom::Sulfur, chemkit::Bond::Double)){
                    setType(index, 6);
                }
            }
            else if(neighborBond->order() == chemkit::Bond::Double){
                if(neighbor->isBondedTo(chemkit::Atom::Nitrogen)){
                    setType(index, 7);
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
                        setType(index, 32, -0.5);
                    }
                    else{
                        setType(index, 7);
                    }
                }
                else if(neighbor->isBondedTo(chemkit::Atom::Carbon) || neighbor->isBondedTo(chemkit::Atom::Hydrogen)){
                    setType(index, 7);
                }
                else{
                    setType(index, 7);
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
                    setType(index, 32);
                }
                else if(neighbor->is(chemkit::Atom::Nitrogen) && oxygenCount == 3 && negativeOxygenCount > 1){
                    setType(index, 32, -1.0/3.0);
                }
                else if(negativeOxygenCount > 1){
                    setType(index, 32, -1.0 / negativeOxygenCount);
                }
                else{
                    setType(index, 32);
                }
            }
            else if(atom->formalCharge() < 0){
                if(neighbor->isBondedTo(chemkit::Atom::Carbon) && neighbor->neighborCount() == 2){
                    setType(index, 35, -1.0);
                }
                else if(neighbor->formalCharge() == 0 && neighbor->neighborCount(chemkit::Atom::Oxygen) == 1){
                    setType(index, 35, -1.0);
                }
                else{
                    setType(index, 32);
                }
            }
            else{
                setType(index, 7);
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
                setType(index, 32); // O-S
            }
            else if(doubleBondedOxygenCount >= 2){
                if(negativeOxygen){
                    setType(index, 32, -1.0/3.0);
                }
                else if(neighbor->valence() == 5 && doubleBondedOxygenCount == 2){
                    setType(index, 32, -0.5);
                }
                else{
                    setType(index, 32); // O2S, O3S, 04S
                }
            }
            else if(singleBondedOxygenCount == 1 && doubleBondedOxygenCount == 1){
                setType(index, 7);
            }
            else if(doubleBondedOxygenCount == 1 &&
                    neighbor->isBondedTo(chemkit::Atom::Sulfur, chemkit::Bond::Double) &&
                    neighbor->valence() == 5){
                setType(index, 32, -0.5); // OSMS
            }
            else{
                setType(index, 7);
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
                    setType(index, 32, -2.0/3.0);
                }
                else{
                    setType(index, 32, -0.5);
                }
            }
            else if(negativeOxygenAndSulfurCount > 1){
                setType(index, 32, -1.0 / negativeOxygenAndSulfurCount);
            }
            else{
                setType(index, 32);
            }
        }
        else if(neighbor->is(chemkit::Atom::Chlorine) && neighbor->neighborCount(chemkit::Atom::Oxygen) == 4){
            setType(index, 32, -0.25); // O4CL
        }
        else if(neighbor->is(chemkit::Atom::Hydrogen) && atom->formalCharge() == -1){
            setType(index, 35, -1.0);
        }
    }

    // two neighbors
    else if(atom->neighborCount() == 2){
        if(atom->neighborCount(chemkit::Atom::Hydrogen) == 2){
            setType(index, 70); // water
        }
        else if(atom->formalCharge() == 1){
            if(atom->isBondedTo(chemkit::Atom::Carbon, chemkit::Bond::Double)){
                setType(index, 51, 1.0);
            }
            else{
                setType(index, 49, 1.0);
            }
        }
        else if(atom->isBondedTo(chemkit::Atom::Nitrogen)){
            setType(index, 6);
        }
        else if(atom->isBondedTo(chemkit::Atom::Sulfur)){
            setType(index, 6);
        }
        else if(atom->isBondedTo(chemkit::Atom::Phosphorus)){
            setType(index, 6);
        }
        else if(atom->isBondedTo(chemkit::Atom::Carbon)){
            setType(index, 6);
        }
        else{
            setType(index, 6);
        }
    }

    // three neighbors
    else if(atom->neighborCount() == 3){
        setType(index, 49, 1.0);
    }
}

void MmffAtomTyper::setSulfurType(int index, const chemkit::Atom *atom)
{
    if(atom->isTerminal()){
        const chemkit::Atom *neighbor = atom->bonds()[0]->otherAtom(atom);

        if(isThiocarboxylate(atom)){
            setType(index, 72, -0.5);
        }
        else if(atom->isBondedTo(chemkit::Atom::Carbon, chemkit::Bond::Double)){
            setType(index, 16);
        }
        else if(atom->isBondedTo(chemkit::Atom::Oxygen, chemkit::Bond::Double)){
            setType(index, 17);
        }
        else if(neighbor->is(chemkit::Atom::Phosphorus)){
            if(neighbor->isBondedTo(chemkit::Atom::Oxygen, chemkit::Bond::Double) && atom->formalCharge() == -1){
                setType(index, 72, -0.5);
            }
            else{
                setType(index, 72); // S-P
            }
        }
        else if(atom->formalCharge() < 0){
            setType(index, 72, -1.0); // SM
        }
        else if(atom->isBondedTo(chemkit::Atom::Sulfur, chemkit::Bond::Double)){
            setType(index, 72, -0.5);
        }
        else{
            setType(index, 72);
        }
    }
    else if(atom->isBondedTo(chemkit::Atom::Nitrogen, chemkit::Bond::Double) &&
            atom->isBondedTo(chemkit::Atom::Oxygen, chemkit::Bond::Double)){
        setType(index, 18);
    }
    else if(atom->isBondedTo(chemkit::Atom::Nitrogen, chemkit::Bond::Double) && atom->neighborCount() == 3){
        setType(index, 17); // >S=N
    }
    else if(atom->isBondedTo(chemkit::Atom::Oxygen, chemkit::Bond::Double) &&
            atom->isBondedTo(chemkit::Atom::Sulfur, chemkit::Bond::Double)){
        setType(index, 73);
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
            setType(index, 17); // S=O
        }
        else if(doubleBondedOxygenCount == 2 && atom->valence() == 5){
            setType(index, 73); // SO2M
        }
        else if(doubleBondedOxygenCount == 1 && atom->isBondedTo(chemkit::Atom::Carbon, chemkit::Bond::Double)){
            setType(index, 74); // =S=O
        }
        else if(doubleBondedOxygenCount >= 2){
            setType(index, 18); // SO2, SO2N, SO3, SO4, =SO2, SNO
        }
        else{
            setType(index, 17); // S=O
        }
    }
    else{
        setType(index, 15);
    }
}

void MmffAtomTyper::setAromaticType(int index, const chemkit::Atom *atom, const chemkit::Ring *ring, int position)
{
    int type = typeNumber(atom);

    // carbon
    if(atom->is(chemkit::Atom::Carbon)){
        if(ring->size() == 5){
            if(type == 57){
                setType(index, 80); // CIM+
            }
            else if(position == 0){
                setType(index, 78); // C5
            }
            else if(position == 1){
                if(type == 64){
                    setType(index, 78); // C5
                }
                else{
                    setType(index, 63); // C5A
                }
            }
            else if(position == 2){
                if(type == 63){
                    setType(index, 78); // C5
                }
                else{
                    setType(index, 64); // C5B
                }
            }
            else{
                setType(index, 78); // C5
            }
        }
        else if(ring->size() == 6){
            setType(index, 37); // CB
        }
    }

    // nitrogen
    else if(atom->is(chemkit::Atom::Nitrogen)){
        if(ring->size() == 5){
            if(type == 62){
                if(ring->atomCount(chemkit::Atom::Nitrogen) == 2){
                    setType(index, 76, -0.5); // N5M
                }
                else if(ring->atomCount(chemkit::Atom::Nitrogen) == 3){
                    setType(index, 76, -1.0/3.0); // N5M
                }
                else if(ring->atomCount(chemkit::Atom::Nitrogen) == 4){
                    setType(index, 76, -1.0/4.0); // N5M
                }
            }
            else if(type == 67){
                setType(index, 82); // N5OX
            }
            else if(type == 54){
                setType(index, 81, 1.0); // N5+
            }
            else if(type == 55){
                setType(index, 81, 0.5); // NIM+
            }
            else if(type == 56){
                setType(index, 81, 1.0/3.0);
            }
            else if(position == 0){
                setType(index, 39); // NPYL
            }
            else if(position == 1){
                if(type == 66){
                    setType(index, 79); // N5
                }
                else{
                    setType(index, 65); // N5A
                }
            }
            else if(position == 2){
                if(type == 65){
                    setType(index, 79); // N5
                }
                else{
                    setType(index, 66); // N5A
                }
            }
            else{
                setType(index, 79); // N5
            }
        }
        else if(ring->size() == 6){
            if(type == 54 || type == 55 || type == 56){
                setType(index, 58, 1.0); // NPD+
            }
            else if(type == 67){
                setType(index, 69); // NPOX
            }
            else{
                if(atom->formalCharge() > 0){
                    setType(index, 58, 1.0); // NPYD+
                }
                else{
                    setType(index, 38); // NPYD
                }
            }
        }
    }

    // oxygen
    else if(atom->is(chemkit::Atom::Oxygen)){
        if(ring->size() == 5){
            setType(index, 59); // OFUR
        }
    }

    // sulfur
    else if(atom->is(chemkit::Atom::Sulfur)){
        if(ring->size() == 5){
            setType(index, 44); // STHI
        }
    }
}
