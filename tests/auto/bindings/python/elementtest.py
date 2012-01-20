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

class ElementTest(unittest.TestCase):
    def test_constructor(self):
        self.assertEqual(chemkit.Element().atomicNumber(), 0)
        self.assertEqual(chemkit.Element(6).atomicNumber(), 6)
        self.assertEqual(chemkit.Element("O").atomicNumber(), 8)
        self.assertEqual(chemkit.Element("").atomicNumber(), 0)

    def test_symbol(self):
        self.assertEqual(chemkit.Element(6).symbol(), "C")
        self.assertEqual(chemkit.Element(2).symbol(), "He")
        self.assertEqual(chemkit.Element().symbol(), "")
        self.assertEqual(chemkit.Element(200).symbol(), "")

    def test_name(self):
        self.assertEqual(chemkit.Element(6).name(), "Carbon")
        self.assertEqual(chemkit.Element(2).name(), "Helium")
        self.assertEqual(chemkit.Element().name(), "")
        self.assertEqual(chemkit.Element(200).name(), "")

    def test_isValid(self):
        self.assertEqual(chemkit.Element().isValid(), False)
        self.assertEqual(chemkit.Element(6).isValid(), True)
        self.assertEqual(chemkit.Element(200).isValid(), False)

if __name__ == '__main__':
    unittest.main()
