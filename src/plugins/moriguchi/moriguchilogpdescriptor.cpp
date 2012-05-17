/******************************************************************************
**
** Copyright (C) 2009-2012 Kyle Lutz <kyle.r.lutz@gmail.com>
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

#include "moriguchilogpdescriptor.h"

#include <chemkit/atom.h>
#include <chemkit/bond.h>
#include <chemkit/ring.h>
#include <chemkit/foreach.h>
#include <chemkit/molecule.h>

// The MoriguchiLogPDescriptor class calculates the Moriguchi logP
// descriptor value for a given molecule.
//
// Reference: [Moriguchi 1992]
MoriguchiLogPDescriptor::MoriguchiLogPDescriptor()
    : chemkit::MolecularDescriptor("moriguchi-logp")
{
    setDimensionality(1);
}

// Returns the Moriguchi logP value for the molecule.
chemkit::Variant MoriguchiLogPDescriptor::value(const chemkit::Molecule *molecule) const
{
    // CX - summation of numbers of carbon and halogen atoms
    //      weighted by C: 1.0, F: 0.5, Cl: 1.0, Br: 1.5, I: 2.0
    chemkit::Real cx = 0;
    foreach(const chemkit::Atom *atom, molecule->atoms()){
        if(atom->is(chemkit::Atom::Carbon)){
            cx += 1.0;
        }
        else if(atom->is(chemkit::Atom::Fluorine)){
            cx += 0.5;
        }
        else if(atom->is(chemkit::Atom::Chlorine)){
            cx += 1.0;
        }
        else if(atom->is(chemkit::Atom::Bromine)){
            cx += 1.5;
        }
        else if(atom->is(chemkit::Atom::Iodine)){
            cx += 2.0;
        }
    }

    // NO - total number of nitrogen and oxygen atoms
    chemkit::Real no = 0;
    foreach(const chemkit::Atom *atom, molecule->atoms()){
        if(atom->is(chemkit::Atom::Nitrogen) ||
           atom->is(chemkit::Atom::Oxygen)){
            no += 1.0;
        }
    }

    // PRX - proximity effect of N/O, X-Y: 2.0, X-A-Y: 1.0
    //       (X, Y: N/O, A: C, S or P) with a correction of
    //       -1 for caroxamide/sulfonamide
    chemkit::Real prx = 0;
    foreach(const chemkit::Atom *atom, molecule->atoms()){
        if(atom->is(chemkit::Atom::Carbon) ||
           atom->is(chemkit::Atom::Sulfur) ||
           atom->is(chemkit::Atom::Phosphorus)){
            if(atom->isBondedTo(chemkit::Atom::Oxygen, chemkit::Bond::Double) &&
               atom->isBondedTo(chemkit::Atom::Oxygen, chemkit::Bond::Single)){
                prx += 2.0;
            }
            else if(atom->isBondedTo(chemkit::Atom::Nitrogen, chemkit::Bond::Double) &&
                    atom->isBondedTo(chemkit::Atom::Nitrogen, chemkit::Bond::Single)){
                prx += 2 * (atom->neighborCount(chemkit::Atom::Nitrogen) - 1);
            }
            else if(atom->isBondedTo(chemkit::Atom::Oxygen, chemkit::Bond::Single) &&
                    atom->isBondedTo(chemkit::Atom::Nitrogen, chemkit::Bond::Single)){
                prx += 2.0;
            }
            else if(atom->isBondedTo(chemkit::Atom::Oxygen,
                                     chemkit::Bond::Double)){
                prx += atom->neighborCount(chemkit::Atom::Nitrogen);
            }
            else if(atom->isBondedTo(chemkit::Atom::Nitrogen, chemkit::Bond::Triple) &&
                    atom->isBondedTo(chemkit::Atom::Nitrogen, chemkit::Bond::Single)){
                prx += 1.0;
            }
            else if(atom->isBondedTo(chemkit::Atom::Oxygen) &&
                    atom->isBondedTo(chemkit::Atom::Nitrogen)){
                prx += 1.0;
            }
        }
        else if(atom->is(chemkit::Atom::Oxygen) &&
                atom->isBondedTo(chemkit::Atom::Nitrogen)){
            prx += 2.0;
        }
    }

    // UB - total number of unsaturated bonds except those in NO2
    chemkit::Real ub = 0;
    foreach(const chemkit::Bond *bond, molecule->bonds()){
        if(!bond->is(chemkit::Bond::Single)){
            ub += 1.0;
        }
    }

    // HB - dummy variable for the presence of intramolecular
    //      hydrogen bond
    chemkit::Real hb = 0;

    // POL - number of aromatic polar substituents (aromatic
    //       substituents excluding Ar-CX2 and Ar-CX=C<, with
    //       X: C or H)
    chemkit::Real pol = 0;
    foreach(const chemkit::Ring *ring, molecule->rings()){
        if(ring->isAromatic()){
            foreach(const chemkit::Bond *exocyclicBond, ring->exocyclicBonds()){
                const chemkit::Atom *substituent = 0;
                if(ring->contains(exocyclicBond->atom1())){
                    substituent = exocyclicBond->atom2();
                }
                else{
                    substituent = exocyclicBond->atom1();
                }

                if(substituent->is(chemkit::Atom::Nitrogen) ||
                   substituent->is(chemkit::Atom::Oxygen) ||
                   substituent->is(chemkit::Atom::Sulfur) ||
                   substituent->is(chemkit::Atom::Fluorine) ||
                   substituent->is(chemkit::Atom::Chlorine) ||
                   substituent->is(chemkit::Atom::Bromine) ||
                   substituent->is(chemkit::Atom::Iodine)){
                    pol += 1.0;
                }
                else if(substituent->is(chemkit::Atom::Carbon)){
                    if(substituent->isBondedTo(chemkit::Atom::Nitrogen,
                                               chemkit::Bond::Double)){
                        pol += 1.0;
                    }
                    else if(substituent->isBondedTo(chemkit::Atom::Oxygen,
                                                    chemkit::Bond::Double)){
                        pol += 1.0;
                    }
                    else if(substituent->isBondedTo(chemkit::Atom::Fluorine) ||
                            substituent->isBondedTo(chemkit::Atom::Chlorine) ||
                            substituent->isBondedTo(chemkit::Atom::Bromine) ||
                            substituent->isBondedTo(chemkit::Atom::Iodine)){
                        pol += 1.0;
                    }
                }
            }
        }
    }

    // AMP - amphoteric property. alpha-aminoacid: 1.0, aminobenzoic
    //       acid: 0.5, pyridinecarboxylic acid: 0.5
    chemkit::Real amp = 0;

    // ALK - dummy variable for alkane, alkene, cycloalkane, or
    //       cycloalkene (hydrocarbons with 0 or 1 double bond)
    chemkit::Real alk = 0;

    // RNG - dummy variable for the presence of ring structures
    //       except benzene and its condensed rings (aromatic,
    //       heteroaromatic, and hydrocarbon rings)
    chemkit::Real rng = 0;
    foreach(const chemkit::Ring *ring, molecule->rings()){
        bool unsaturated = false;

        foreach(const chemkit::Bond *bond, ring->bonds()){
            if(bond->is(chemkit::Bond::Double)){
                unsaturated = true;
            }
        }

        if(!unsaturated || ring->isHeterocycle() || ring->size() != 6){
            rng = 1.0;
            break;
        }
    }

    // QN - quaternary nitrogen: >N+<, 1.0; N oxide, 0.5
    chemkit::Real qn = 0;
    foreach(const chemkit::Atom *atom, molecule->atoms()){
        if(atom->is(chemkit::Atom::Nitrogen) &&
           atom->neighborCount() == 4){
            qn = 1.0;
            break;
        }
    }

    // NO2 - number of nitro groups
    chemkit::Real no2 = 0;
    foreach(const chemkit::Atom *atom, molecule->atoms()){
        if(atom->is(chemkit::Atom::Nitrogen) &&
           atom->neighborCount() == 3 &&
           atom->isBondedTo(chemkit::Atom::Oxygen, chemkit::Bond::Single) &&
           atom->isBondedTo(chemkit::Atom::Oxygen, chemkit::Bond::Double)){
            no2 += 1.0;

            // remove one from unsaturated bond count
            ub -= 1.0;
        }
    }

    // NCS - isothiocyanato (-N=C=S), 1.0;
    //       thiocyanato (-S-C#N), 0.5
    chemkit::Real ncs = 0;
    foreach(const chemkit::Atom *atom, molecule->atoms()){
        if(atom->is(chemkit::Atom::Carbon)){
            if(atom->isBondedTo(chemkit::Atom::Sulfur, chemkit::Bond::Double) &&
               atom->isBondedTo(chemkit::Atom::Nitrogen, chemkit::Bond::Single)){
                ncs += 1.0;
            }
            else if(atom->isBondedTo(chemkit::Atom::Sulfur, chemkit::Bond::Single) &&
                    atom->isBondedTo(chemkit::Atom::Nitrogen, chemkit::Bond::Triple)){
                ncs += 0.5;
            }
        }
    }

    // BLM - dummy variable for the presence of beta-lactam
    chemkit::Real blm = 0;
    foreach(const chemkit::Ring *ring, molecule->rings()){
        if(ring->size() == 4){
            bool containsNitrogen = false;
            bool containsCarbonyl = false;

            foreach(const chemkit::Atom *atom, ring->atoms()){
                if(atom->is(chemkit::Atom::Nitrogen)){
                    containsNitrogen = true;
                }
                else if(atom->is(chemkit::Atom::Carbon)){
                    if(atom->isBondedTo(chemkit::Atom::Oxygen,
                                        chemkit::Bond::Double) &&
                       atom->isBondedTo(chemkit::Atom::Nitrogen,
                                        chemkit::Bond::Single)){
                        containsCarbonyl = true;
                    }
                }
            }

            if(containsNitrogen && containsCarbonyl){
                blm = 1.0;
            }
        }
    }

    // equation 4 in [Moriguchi 1992]
    return 1.244 * pow(cx, 0.6) -
           1.017 * pow(no, 0.9) +
           0.406 * prx -
           0.145 * pow(ub, 0.8) +
           0.511 * hb +
           0.268 * pol -
           2.215 * amp +
           0.912 * alk -
           0.392 * rng -
           3.684 * qn +
           0.474 * no2 +
           1.582 * ncs +
           0.773 * blm -
           1.041;
}
