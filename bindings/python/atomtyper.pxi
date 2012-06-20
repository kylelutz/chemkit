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

from string cimport string

from atomtyper cimport _AtomTyper
from atomtyper cimport create as _AtomTyper_create
from atomtyper cimport typers as _AtomTyper_typers
from atomtyper cimport assignAtomTypes as _AtomTyper_assignAtomTypes
from atomtyper cimport isCarbonylCarbon as _AtomTyper_isCarbonylCarbon
from atomtyper cimport isCarbonylOxygen as _AtomTyper_CarbonylOxygen
from atomtyper cimport isHalogen as _AtomTyper_Halogen
from atomtyper cimport isHydrogenDonor as _AtomTyper_HydrogenDonor
from atomtyper cimport isHydrogenAcceptor as _AtomTyper_HydrogenAcceptor
from atomtyper cimport isHydroxylHydrogen as _AtomTyper_HydroxylHydrogen
from atomtyper cimport isHydroxylOxygen as _AtomTyper_HydroxylOxygen
from atomtyper cimport isNitrileCarbon as _AtomTyper_NitrileCarbon
from atomtyper cimport isNitrileNitrogen as _AtomTyper_NitrileNitrogen
from atomtyper cimport isNitroOxygen as _AtomTyper_NitroOxygen
from atomtyper cimport isNitroNitrogen as _AtomTyper_NitroNitrogen
from atomtyper cimport isPolarHydrogen as _AtomTyper_PolarHydrogen
from atomtyper cimport isNonpolarHydrogen as _AtomTyper_NonpolarHydrogen
from atomtyper cimport isThiolHydrogen as _AtomTyper_ThiolHydrogen
from atomtyper cimport isThiolSulfur as _AtomTyper_ThiolSulfur

cdef class AtomTyper:
    """The AtomTyper class assigns symbolic types to atoms."""

    cdef _AtomTyper *_typer

    ### Construction and Destruction ##########################################
    def __cinit__(self):
        self._typer = NULL

    def __init__(self):
        raise TypeError("AtomTyper objects are not constructible, "
                        "use the AtomTyper.create() method")

    def __dealloc__(self):
        del self._typer

    ### Properties ############################################################
    def name(self):
        """Returns the name of the atom typer."""

        cdef string name = self._typer.name()

        return name.c_str()

    def setMolecule(self, Molecule molecule):
        """Sets the molecule."""

        self._typer.setMolecule(molecule._molecule)

    def molecule(self):
        """Returns the molecule."""

        return Molecule_fromPointer(self._typer.molecule())

    ### Types #################################################################
    def type(self, Atom atom):
        """Returns the type for the atom."""

        cdef string type = self._typer.type(atom._atom)

        return type.c_str()

    ### Predicates ############################################################
    @classmethod
    def isCarbonylCarbon(cls, Atom atom):
        """Returns True if the atom is the Carbon in a carbonyl group."""

        return _AtomTyper_isCarbonylCarbon(atom._atom)

    @classmethod
    def isCarbonylOxygen(cls, Atom atom):
        """Returns True if the atom is the Oxygen in a carbonyl group."""

        return _AtomTyper_CarbonylOxygen(atom._atom)

    @classmethod
    def isHalogen(cls, Atom atom):
        """Returns True if the atom is a halogen."""

        return _AtomTyper_Halogen(atom._atom)

    @classmethod
    def isHydrogenDonor(cls, Atom atom):
        """Returns True if the atom is a hydrogen donor."""

        return _AtomTyper_HydrogenDonor(atom._atom)

    @classmethod
    def isHydrogenAcceptor(cls, Atom atom):
        """Returns True if the atom is a hydrogen acceptor."""

        return _AtomTyper_HydrogenAcceptor(atom._atom)

    @classmethod
    def isHydroxylHydrogen(cls, Atom atom):
        """Returns True if the atom is the Hydrogen in a hydroxyl group."""

        return _AtomTyper_HydroxylHydrogen(atom._atom)

    @classmethod
    def isHydroxylOxygen(cls, Atom atom):
        """Returns True if the atom is the Oxygen in a hydroxyl group."""

        return _AtomTyper_HydroxylOxygen(atom._atom)

    @classmethod
    def isNitrileCarbon(cls, Atom atom):
        """Returns True if the atom is the Carbon in a nitrile group."""

        return _AtomTyper_NitrileCarbon(atom._atom)

    @classmethod
    def isNitrileNitrogen(cls, Atom atom):
        """Returns True if the atom is the Nitrogen in a nitrile group."""

        return _AtomTyper_NitrileNitrogen(atom._atom)

    @classmethod
    def isNitroOxygen(cls, Atom atom):
        """Returns True if the atom is an Oxygen in a nitro group."""

        return _AtomTyper_NitroOxygen(atom._atom)

    @classmethod
    def isNitroNitrogen(cls, Atom atom):
        """Returns True if the atom is the Nitrogen in a nitro group."""

        return _AtomTyper_NitroNitrogen(atom._atom)

    @classmethod
    def isPolarHydrogen(cls, Atom atom):
        """Returns True if the atom is a polar Hydrogen."""

        return _AtomTyper_PolarHydrogen(atom._atom)

    @classmethod
    def isNonpolarHydrogen(cls, Atom atom):
        """Returns True if the atom is a non-polar Hydrogen."""

        return _AtomTyper_NonpolarHydrogen(atom._atom)

    @classmethod
    def isThiolHydrogen(cls, Atom atom):
        """Returns True if the atom is the Hydrogen in a thiol group."""

        return _AtomTyper_ThiolHydrogen(atom._atom)

    @classmethod
    def isThiolSulfur(cls, Atom atom):
        """Returns True if the atom is the Sulfur in a thiol group."""

        return _AtomTyper_ThiolSulfur(atom._atom)

    ### Static Methods ########################################################
    @classmethod
    def create(cls, char *name):
        """Creates a new atom typer."""

        cdef _AtomTyper *_typer = _AtomTyper_create(name)
        if _typer is NULL:
            return None

        cdef AtomTyper typer = AtomTyper.__new__(AtomTyper)
        typer._typer = _typer
        return typer

    @classmethod
    def typers(cls):
        """Returns a list of the names of the supported atom typers."""

        cdef vector[string] _typers = _AtomTyper_typers()

        typers = []
        for i in range(_typers.size()):
            typers.append(_typers[i].c_str())

        return typers

    @classmethod
    def assignAtomTypes(cls, Molecule molecule, char *typer):
        """Assigns types to the atoms in the molecule using typer."""

        return _AtomTyper_assignAtomTypes(molecule._molecule, <string>(typer))
