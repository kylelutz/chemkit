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

from string cimport string

from atom cimport _Atom

cdef class Point3
cdef class Molecule

cdef class Atom:
    """The Atom class represents an atom in a molecule.

    Atom objects are created with the
    :func:`~chemkit.Molecule.addAtom` method and destroyed with the
    :func:`~chemkit.Molecule.removeAtom` method.

    """

    cdef _Atom *_atom

    ### Construction and Destruction ##########################################
    def __cinit__(self):
        self._atom = NULL

    def __init__(self):
        raise TypeError("Atom objects are not constructible")

    ### Properties ############################################################
    def element(self):
        """Returns the element for the atom."""

        return Element(self.atomicNumber())

    def setAtomicNumber(self, int atomicNumber):
        """Sets the atomic number for the atom."""

        self._atom.setAtomicNumber(atomicNumber)

    def atomicNumber(self):
        """Returns the atomic number of the atom."""

        return self._atom.atomicNumber()

    def symbol(self):
        """Returns the element symbol for the atom."""

        return self._atom.symbol().c_str()

    def name(self):
        """Returns the name of the element for the atom."""

        return self._atom.name().c_str()

    def setType(self, char *type):
        """Sets the symbolic type for the atom."""

        self._atom.setType(<string>(type))

    def type(self):
        """Returns the symbolic type for the atom."""

        cdef string type = self._atom.type()

        return type.c_str()

    def formalCharge(self):
        """Returns the formal charge of the atom."""

        return self._atom.formalCharge()

    def setPartialCharge(self, double charge):
        """Sets the partial charge for the atom."""

        self._atom.setPartialCharge(charge)

    def partialCharge(self):
        """Returns the partial charge of the atom."""

        return self._atom.partialCharge()

    def mass(self):
        """Returns the mass of the atom."""

        return self._atom.mass()

    def electronegativity(self):
        """Returns the electronegativity of the atom."""

        return self._atom.electronegativity()

    def covalentRadius(self):
        """Returns the covalent radius of the atom."""

        return self._atom.covalentRadius()

    def vanDerWaalsRadius(self):
        """Returns the Van der Waals radius of the atom."""

        return self._atom.vanDerWaalsRadius()

    def fragment(self):
        """Returns the fragment that the atom is a part of."""

        return Fragment_fromPointer(self._atom.fragment())

    def molecule(self):
        """Returns the molecule that the atom is a part of."""

        return Molecule_fromPointer(self._atom.molecule())

    def index(self):
        """Returns the index of the atom in the molecule."""

        return self._atom.index()

    ### Structure #############################################################
    def bond(self, index):
        """Returns the bond at index for the atom."""

        return Bond_fromPointer(self._atom.bond(index))

    def bonds(self):
        """Returns a list of bonds that the atom is a part of."""

        bonds = []
        for i in range(self.bondCount()):
            bonds.append(self.bond(i))

        return bonds

    def bondCount(self):
        """Returns the number of bonds that the atom is a part of."""

        return self._atom.bondCount()

    def valence(self):
        """Returns the valence of the atom."""

        return self._atom.valence()

    def bondTo(self, Atom neighbor):
        """Returns the bond from the atom to neighbor."""

        return Bond_fromPointer(self._atom.bondTo(neighbor._atom))

    def neighbor(self, int index):
        """Returns the neighbor at index."""

        return Atom_fromPointer(self._atom.neighbor(index))

    def neighbors(self):
        """Returns a list of neighboring atoms."""

        neighbors = []
        for i in range(self.neighborCount()):
            neighbors.append(self.neighbor(i))

        return neighbors

    def neighborCount(self):
        """Returns the number of neighboring atoms."""

        return self._atom.neighborCount()

    def isBondedTo(self, Atom atom):
        """Returns True if the atom is bonded to atom."""

        return self._atom.isBondedTo(atom._atom)

    def isConnectedTo(self, Atom atom):
        """Returns True if the atom is connected to atom."""

        return self._atom.isConnectedTo(atom._atom)

    def isTerminal(self):
        """Returns True if the atom is terminal."""

        return self._atom.isTerminal()

    def isTerminalHydrogen(self):
        """Returns True if the atom is terminal and is a hydrogen."""

        return self._atom.isTerminalHydrogen()

    ### Ring Perception #######################################################
    def ring(self, index):
        """Returns the ring at index for the atom."""

        return Ring_fromPointer(self._atom.ring(index))

    def rings(self):
        """Returns a list of rings that the atom is a part of."""

        rings = []
        for i in range(self.ringCount()):
            rings.append(self.ring(i))

        return rings

    def ringCount(self):
        """Returns the number of rings that the atom is a part of."""

        return self._atom.ringCount()

    def isInRing(self, int size = 0):
        """Returns True if the atom is in a ring."""

        if size != 0:
            return self._atom.isInRing(size)
        else:
            return self._atom.isInRing()

    def smallestRing(self):
        """Returns the smallest ring that the atom is a part of."""

        return Ring_fromPointer(self._atom.smallestRing())

    def isAromatic(self):
        """Returns True if the atom is in an aromatic ring."""

        return self._atom.isAromatic()

    ### Geometry ##############################################################
    def setPosition(self, x, y, z):
        """Sets the position of the atom to (x, y, z)."""

        self._atom.setPosition(x, y, z)

    def position(self):
        """Returns the position of the atom."""

        cdef _Point3 p = self._atom.position()

        return Point3(p.x(), p.y(), p.z())

    def x(self):
        """Returns the x coordinate of the atoms position."""

        return self._atom.x()

    def y(self):
        """Returns the y coordinate of the atoms position."""

        return self._atom.y()

    def z(self):
        """Returns the z coordinate of the atoms position."""

        return self._atom.z()

    def distance(self, Atom atom):
        """Returns the distance from the atom to atom."""

        return self._atom.distance(atom._atom)

    ### Operators #############################################################
    def __richcmp__(Atom self, Atom atom, int op):
        if op == 2:
            return self._atom == atom._atom
        else:
            return False

    # enum AtomName
    Hydrogen = 1
    Helium = 2
    Lithium = 3
    Beryllium = 4
    Boron = 5
    Carbon = 6
    Nitrogen = 7
    Oxygen = 8
    Fluorine = 9
    Neon = 10
    Sodium = 11
    Magnesium = 12
    Aluminum = 13
    Silicon = 14
    Phosphorus = 15
    Sulfur = 16
    Chlorine = 17
    Argon = 18
    Potassium = 19
    Calcium = 20
    Scandium = 21
    Titanium = 22
    Vanadium = 23
    Chromium = 24
    Manganese = 25
    Iron = 26
    Cobalt = 27
    Nickel = 28
    Copper = 29
    Zinc = 30
    Gallium = 31
    Germanium = 32
    Arsenic = 33
    Selenium = 34
    Bromine = 35
    Krypton = 36
    Rubidium = 37
    Strontium = 38
    Yttrium = 39
    Zirconium = 40
    Niobium = 41
    Molybdenum = 42
    Technetium = 43
    Ruthenium = 44
    Rhodium = 45
    Palladium = 46
    Silver = 47
    Cadmium = 48
    Indium = 49
    Tin = 50
    Antimony = 51
    Tellurium = 52
    Iodine = 53
    Xenon = 54
    Cesium = 55
    Barium = 56
    Lanthanum = 57
    Cerium = 58
    Praseodymium = 59
    Neodymium = 60
    Promethium = 61
    Samarium = 62
    Europium = 63
    Gadolinium = 64
    Terbium = 65
    Dysprosium = 66
    Holmium = 67
    Erbium = 68
    Thulium = 69
    Ytterbium = 70
    Lutetium = 71
    Hafnium = 72
    Tantalum = 73
    Tungsten = 74
    Rhenium = 75
    Osmium = 76
    Iridium = 77
    Platinum = 78
    Gold = 79
    Mercury = 80
    Thallium = 81
    Lead = 82
    Bismuth = 83
    Polonium = 84
    Astatine = 85
    Radon = 86
    Francium = 87
    Radium = 88
    Actinium = 89
    Thorium = 90
    Protactinium = 91
    Uranium = 92
    Neptunium = 93
    Plutonium = 94
    Americium = 95
    Curium = 96
    Berkelium = 97
    Californium = 98
    Einsteinium = 99
    Fermium = 100
    Mendelevium = 101
    Nobelium = 102
    Lawrencium = 103
    Rutherfordium = 104
    Dubnium = 105
    Seaborgium = 106
    Bohrium = 107
    Hassium = 108
    Meitnerium = 109

cdef Atom Atom_fromPointer(_Atom *_atom):
    cdef Atom atom = Atom.__new__(Atom)
    atom._atom = _atom
    return atom

