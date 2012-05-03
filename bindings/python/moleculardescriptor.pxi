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
from libcpp.vector cimport vector

from moleculardescriptor cimport _MolecularDescriptor
from moleculardescriptor cimport create as _MolecularDescriptor_create
from moleculardescriptor cimport descriptors as _MolecularDescriptor_descriptors

cdef class MolecularDescriptor:
    cdef _MolecularDescriptor *_descriptor

    ### Construction and Destruction ##########################################
    def __cinit__(self):
        self._descriptor = NULL

    def __init__(self):
        raise TypeError("MolecularDescriptor objects are not constructible, use the MolecularDescriptor.create() method")

    def __dealloc__(self):
        del self._descriptor

    ### Properties ############################################################
    def name(self):
        """Returns the name of the descriptor."""

        return self._descriptor.name().c_str()

    def dimensionality(self):
        """Returns the dimensionality of the descriptor."""

        return self._descriptor.dimensionality()

    ### Descriptor ############################################################
    def value(self, Molecule molecule):
        """Returns the descriptor value for the molecule."""

        cdef _Variant value = self._descriptor.value(molecule._molecule)

        return value.toDouble()

    ### Static Methods ########################################################
    @classmethod
    def create(cls, char *name):
        """Creates a new molecular descriptor."""

        cdef _MolecularDescriptor *_descriptor = _MolecularDescriptor_create(name)
        if _descriptor is NULL:
            return None

        cdef MolecularDescriptor descriptor = MolecularDescriptor.__new__(MolecularDescriptor)
        descriptor._descriptor = _descriptor

        return descriptor

    @classmethod
    def descriptors(cls):
        """Returns a list of the names of the supported molecular descriptors."""

        cdef vector[string] _descriptors = _MolecularDescriptor_descriptors()

        descriptors = []
        for i in range(_descriptors.size()):
            descriptors.append(_descriptors[i].c_str())

        return descriptors

