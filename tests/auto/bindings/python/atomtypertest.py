###############################################################################
##
## Copyright (C) 2012 Kitware, Inc.
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

class AtomTyperTest(unittest.TestCase):
    def test_name(self):
        typer = chemkit.AtomTyper.create("uff")
        self.assertEqual(typer.name(), "uff")

    def test_benzene_mmff(self):
        molecule = chemkit.Molecule("c1ccccc1", "smiles")
        self.assertEqual(molecule.formula(), "C6H6")

        typer = chemkit.AtomTyper.create("mmff")
        typer.setMolecule(molecule)

        for atom in molecule.atoms():
            if atom.symbol() == "C":
                self.assertEqual(typer.type(atom), "37")
            elif atom.symbol() == "H":
                self.assertEqual(typer.type(atom), "5")

    def test_assign_types(self):
        molecule = chemkit.Molecule("CCO", "smiles")
        self.assertEqual(molecule.formula(), "C2H6O")

        for atom in molecule.atoms():
            self.assertEqual(atom.type(), "")

        chemkit.AtomTyper.assignAtomTypes(molecule, "element-name")

        for atom in molecule.atoms():
            if atom.symbol() == "C":
                self.assertEqual(atom.type(), "Carbon")
            elif atom.symbol() == "H":
                self.assertEqual(atom.type(), "Hydrogen")
            elif atom.symbol() == "O":
                self.assertEqual(atom.type(), "Oxygen")

if __name__ == '__main__':
    unittest.main()
