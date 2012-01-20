###############################################################################
##
## Copyright (C) 2009-2012 Kyle Lutz <kyle.r.lutz@gmail.com>
## All rights reserved.
##
## This file is a part of the chemkit project. For more information
## see <http://www.chemkit.org>.
##
## Redistribution and use in source and binary forms, with or without
## modification, are permitted provided that the following conditions
## are met:
##
##   * Redistributions of source code must retain the above copyright
##     notice, this list of conditions and the following disclaimer.
##   * Redistributions in binary form must reproduce the above copyright
##     notice, this list of conditions and the following disclaimer in the
##     documentation and/or other materials provided with the distribution.
##   * Neither the name of the chemkit project nor the names of its
##     contributors may be used to endorse or promote products derived
##     from this software without specific prior written permission.
##
## THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
## "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
## LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
## A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
## OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
## SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
## LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
## DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
## THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
## (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
## OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
##
###############################################################################

import chemkit
import unittest

class MoleculeTest(unittest.TestCase):
    def test_atom(self):
        molecule = chemkit.Molecule()
        C1 = molecule.addAtom("C")
        self.assertEqual(C1.atomicNumber(), 6)
        C2 = molecule.addAtom(6)
        self.assertEqual(C2.atomicNumber(), 6)
        C3 = molecule.addAtom(chemkit.Atom.Carbon)
        self.assertEqual(C3.atomicNumber(), 6)
        self.assertEqual(molecule.atomCount(), 3)
        self.assertEqual(molecule.formula(), "C3")

    def test_bond(self):
        molecule = chemkit.Molecule()
        self.assertEqual(molecule.bondCount(), 0)
        C1 = molecule.addAtom("C")
        C2 = molecule.addAtom("C")
        C3 = molecule.addAtom("C")
        C1_C2 = molecule.addBond(C1, C2)
        self.assertEqual(molecule.bondCount(), 1)
        C2_C3 = molecule.addBond(C2, C3, 2)
        self.assertEqual(molecule.bondCount(), 2)
        self.assertEqual(C2_C3.order(), 2)

    def test_formula(self):
        # empty molecule
        empty = chemkit.Molecule()
        self.assertEqual(empty.formula(), "")

        # ethanol molecule
        ethanol = chemkit.Molecule()
        C1 = ethanol.addAtom("C")
        C2 = ethanol.addAtom("C")
        O3 = ethanol.addAtom("O")
        ethanol.addBond(C1, C2)
        ethanol.addBond(C2, O3)

        for atom in ethanol.atoms():
            while(atom.formalCharge() < 0):
                hydrogen = ethanol.addAtom("H")
                ethanol.addBond(atom, hydrogen)

        self.assertEqual(ethanol.formula(), "C2H6O")

    def test_rings(self):
        # benzene molecule
        benzene = chemkit.Molecule("c1ccccc1", "smiles")
        self.assertEqual(benzene.formula(), "C6H6")
        self.assertEqual(benzene.ringCount(), 1)
        ring = benzene.ring(0)
        self.assertEqual(ring.size(), 6)

if __name__ == '__main__':
    unittest.main()
