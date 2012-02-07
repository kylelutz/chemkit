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

from bond cimport _Bond

cdef class Bond:
    """The Bond class represents a bond between two atoms.

    Bond objects are created with the
    :func:`~chemkit.Molecule.addBond` method and destroyed with the
    :func:`~chemkit.Molecule.removeBond` method.

    """

    cdef _Bond *_bond

    # enum BondType
    Single = 1
    Double = 2
    Triple = 3
    Quadruple = 4

    ### Construction and Destruction ##########################################
    def __cinit__(self):
        self._bond = NULL

    def __init__(self):
        raise TypeError("Bond objects are not constructible")

    ### Properties ############################################################
    def atom(self, int index):
        """Returns the atom atom index in the bond."""

        return Atom_fromPointer(self._bond.atom(index))

    def atom1(self):
        """Returns the first atom in the bond."""

        return Atom_fromPointer(self._bond.atom1())

    def atom2(self):
        """Returns the second atom in the bond."""

        return Atom_fromPointer(self._bond.atom2())

    def atoms(self):
        """Returns a list containing both atoms in the bond."""

        return [self.atom1(), self.atom2()]

    def setOrder(self, int order):
        """Sets the bond order for the bond to order."""

        self._bond.setOrder(order)

    def order(self):
        """Returns the bond order of the bond."""

        return self._bond.order()

    def polarity(self):
        """Returns the polarity of the bond."""

        return self._bond.polarity()

    def index(self):
        """Returns the index of the bond."""

        return self._bond.index()

    def fragment(self):
        """Returns the fragment that the bond is a part of."""

        return Fragment_fromPointer(self._bond.fragment())

    def molecule(self):
        """Returns the molecule that the bond is a part of."""

        return Molecule_fromPointer(self._bond.molecule())

    ### Ring Perception #######################################################
    def ring(self, index):
        """Returns the ring at index for the bond."""

        return Ring_fromPointer(self._bond.ring(index))

    def rings(self):
        """Returns a list of rings that the bond is a part of."""

        rings = []
        for i in range(self.ringCount()):
            rings.append(self.ring(i))

        return rings

    def ringCount(self):
        """Returns the number of rings that the bond is a part of."""

        return self._bond.ringCount()

    def isInRing(self, int size = 0):
        """Returns True if the bond is in a ring."""

        if size != 0:
            return self._bond.isInRing(size)
        else:
            return self._bond.isInRing()

    def smallestRing(self):
        """Returns the smallest ring that the bond is a part of."""

        return Ring_fromPointer(self._bond.smallestRing())

    def isAromatic(self):
        """Returns True if the bond is in an aromatic ring."""

        return self._bond.isAromatic()

    ### Geometry ##############################################################
    def length(self):
        """Returns the length of the bond in Angstroms."""

        return self._bond.length()

cdef Bond Bond_fromPointer(_Bond *_bond):
    cdef Bond bond = Bond.__new__(Bond)
    bond._bond = _bond
    return bond

