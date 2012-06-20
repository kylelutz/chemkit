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

dataPath = "../../../data/"

class MoleculeFileTest(unittest.TestCase):
    def test_fileName(self):
        f = chemkit.MoleculeFile()
        self.assertEqual(f.fileName(), "")

        f.setFileName("methane.sdf")
        self.assertEqual(f.fileName(), "methane.sdf")

        f.setFileName("")
        self.assertEqual(f.fileName(), "")

        f = chemkit.MoleculeFile("ethanol.cml")
        self.assertEqual(f.fileName(), "ethanol.cml")

    def test_methanol(self):
        f = chemkit.MoleculeFile(dataPath + "methanol.sdf")
        ok = f.read()
        if not ok:
            print f.errorString()
        self.assertTrue(ok)
        self.assertEqual(f.moleculeCount(), 1)

        molecule = f.molecule(0)
        self.assertIsNotNone(molecule)
        self.assertEqual(molecule.formula(), "CH4O")

    def test_readString(self):
        f = chemkit.MoleculeFile()
        f.setFormat("xyz")

        # string containing the file's data
        data = "5\n"
        data += "methane\n"
        data += "C 	 0.000000 	 0.000000 	 0.000000\n"
        data += "H 	 0.000000 	 0.000000 	 1.089000\n"
        data += "H 	 1.026719 	 0.000000 	-0.363000\n"
        data += "H 	-0.513360 	-0.889165 	-0.363000\n"
        data += "H 	-0.513360 	 0.889165 	-0.363000\n"

        ok = f.readString(data)
        if not ok:
            print f.errorString()
        self.assertTrue(ok)
        self.assertEqual(f.moleculeCount(), 1)

        molecule = f.molecule(0)
        self.assertIsNotNone(molecule)
        self.assertEqual(molecule.formula(), "CH4")

if __name__ == '__main__':
    unittest.main()
