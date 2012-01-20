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

from libcpp cimport bool
from string cimport string

from element cimport _Element

cdef class Element:
    """The Element class represents an element.

    Element objects can be created using either their atomic number
    or their symbol. Either of the following lines will create a
    Carbon element::

        # create carbon element from its atomic number
        carbon = chemkit.Element(6)

        # create carbon element from its symbol
        carbon = chemkit.Element("C")

    """

    cdef _Element *_element

    ### Construction and Destruction ##########################################
    def __init__(self, symbolOrNumber = None):
        # construct empty element
        if symbolOrNumber is None:
            self._element = new _Element()
            return

        # try constructing from atomic number
        try:
            self._element = new _Element(<int?>(symbolOrNumber))
            return
        except:
            pass

        # try constructing from symbol
        try:
            self._element = new _Element(<char*?>(symbolOrNumber))
            return
        except:
            pass

        self._element = new _Element()

    ### Properties ############################################################
    def setAtomicNumber(self, int atomicNumber):
        """Sets the atomic number for the element."""

        self._element.setAtomicNumber(atomicNumber)

    def atomicNumber(self):
        """Returns the atomic number for the element."""

        return self._element.atomicNumber()

    def symbol(self):
        """Returns the symbol for the element."""

        return self._element.symbol().c_str()

    def name(self):
        """Returns the name for the element."""

        return self._element.name().c_str()

    def period(self):
        """Returns the period in the periodic table that the element
           belongs to."""

        return self._element.period()

    def mass(self):
        """Returns the mass of the element."""

        return self._element.mass()

    def electronegativity(self):
        """Returns the electronegativity of the element."""

        return self._element.electronegativity()

    def covalentRadius(self):
        """Returns the covalent radius of the element."""

        return self._element.covalentRadius()

    def vanDerWaalsRadius(self):
        """Returns the Van der Waals radius of the element."""

        return self._element.vanDerWaalsRadius()

    def expectedValence(self):
        """Returns the expected valence for the element."""

        return self._element.expectedValence()

    def isValid(self):
        """Returns True if the element is valid."""

        return self._element.isValid()

    def isMetal(self):
        """Returns True if the element is a metal."""

        return self._element.isMetal()

    def isNonmetal(self):
        """Returns True if the element is a nonmetal."""

        return self._element.isNonmetal()

