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
from libcpp.vector cimport vector

from fragment cimport _Fragment

cdef class Fragment:
    """The Fragment class represents a group of connected atoms."""

    cdef _Fragment *_fragment

    ### Construction and Destruction ##########################################
    def __cinit__(self):
        self._fragment = NULL

    def __init__(self):
        raise TypeError("Fragment objects are not constructible")

    ### Properties ############################################################
    def size(self):
        """Returns the number of atoms in the fragment."""

        return self._fragment.size()

    def molecule(self):
        """Returns the molecule that the fragment is a part of."""

        return Molecule_fromPointer(self._fragment.molecule())

    ### Structure #############################################################
    def atom(self, int index):
        """Returns the atom at index in the fragment."""

        return Atom_fromPointer(self._fragment.atom(index))

    def atoms(self):
        """Returns a list of the atoms in the fragment."""

        atoms = []
        for i in range(self.atomCount()):
            atoms.append(self.atom(i))

        return atoms

    def atomCount(self):
        """Returns the number of atoms in the fragment."""

        return self._fragment.atomCount()

    def contains(self, Atom atom):
        """Returns True if the fragment contains atom."""

        return self._fragment.contains(atom._atom)

    def bondCount(self):
        """Returns the number of bonds in the fragment."""

        return self._fragment.bondCount()

cdef Fragment Fragment_fromPointer(_Fragment *_fragment):
    cdef Fragment fragment = Fragment.__new__(Fragment)
    fragment._fragment = _fragment
    return fragment

