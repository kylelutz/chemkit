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

from libcpp cimport bool
from string cimport string
from libcpp.vector cimport vector

cdef extern from "chemkit/atom.h" namespace "chemkit":
    cdef cppclass _Atom "chemkit::Atom"

cdef extern from "chemkit/molecule.h" namespace "chemkit":
    cdef cppclass _Molecule "chemkit::Molecule"

cdef extern from "chemkit/atomtyper.h" namespace "chemkit":
    cdef cppclass _AtomTyper "chemkit::AtomTyper":
        # properties
        string name()
        void setMolecule(_Molecule *molecule)
        _Molecule* molecule()

        # types
        string type(_Atom *atom)

# static methods
cdef extern from "chemkit/atomtyper.h" namespace "chemkit::AtomTyper":
    _AtomTyper* create(char *name)
    vector[string] typers()
    bool assignAtomTypes(_Molecule *molecule, string typer)

    # predicates
    bool isCarbonylCarbon(_Atom *atom)
    bool isCarbonylOxygen(_Atom *atom)
    bool isHalogen(_Atom *atom)
    bool isHydrogenDonor(_Atom *atom)
    bool isHydrogenAcceptor(_Atom *atom)
    bool isHydroxylHydrogen(_Atom *atom)
    bool isHydroxylOxygen(_Atom *atom)
    bool isNitrileCarbon(_Atom *atom)
    bool isNitrileNitrogen(_Atom *atom)
    bool isNitroOxygen(_Atom *atom)
    bool isNitroNitrogen(_Atom *atom)
    bool isPolarHydrogen(_Atom *atom)
    bool isNonpolarHydrogen(_Atom *atom)
    bool isThiolHydrogen(_Atom *atom)
    bool isThiolSulfur(_Atom *atom)
