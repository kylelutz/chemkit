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

cdef extern from "chemkit/atom.h" namespace "chemkit":
    cdef cppclass _Atom "chemkit::Atom"

cdef extern from "chemkit/bond.h" namespace "chemkit":
    cdef cppclass _Bond "chemkit::Bond"

cdef extern from "chemkit/element.h" namespace "chemkit":
    cdef cppclass _Element "chemkit::Element"

cdef extern from "chemkit/fragment.h" namespace "chemkit":
    cdef cppclass _Fragment "chemkit::Fragment"

cdef extern from "chemkit/ring.h" namespace "chemkit":
    cdef cppclass _Ring "chemkit::Ring"

cdef extern from "chemkit/variant.h" namespace "chemkit":
    cdef cppclass _Variant "chemkit::Variant"

cdef extern from "chemkit/molecule.h" namespace "chemkit":
    cdef cppclass _Molecule "chemkit::Molecule":
        # construction and destruction
        _Molecule()
        _Molecule(char* formula, char* format)

        # properties
        void setName(char* name)
        string name()
        string formula()
        string formula(char* format)
        _Variant descriptor(char* name)
        int size()
        bool isEmpty()
        double mass()
        _Variant data(char* name)

        # structure
        _Atom* addAtom(_Element element)
        void removeAtom(_Atom *atom)
        _Atom* atom(int index)
        int atomCount()
        bool contains(_Atom *atom)
        bool contains(_Element element)
        _Bond* addBond(_Atom *a, _Atom *b, int order)
        void removeBond(_Bond *bond)
        _Bond* bond(int index)
        int bondCount()

        # ring perception
        _Ring* ring(int index)
        int ringCount()

        # fragment perception
        _Fragment* fragment(int index)
        int fragmentCount()
        bool isFragmented()
        void removeFragment(_Fragment *fragment)

