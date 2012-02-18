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

from forcefield cimport _ForceField
from forcefield cimport create as _ForceField_create
from forcefield cimport forceFields as _ForceField_forceFields

cdef class ForceField:
    cdef _ForceField *_forceField

    ### Construction and Destruction ##########################################
    def __cinit__(self):
        self._forceField = NULL

    def __init__(self):
        raise TypeError("ForceField objects are not constructible, use the ForceField.create() method")

    def __dealloc__(self):
        del self._forceField

    ### Properties ############################################################
    def name(self):
        """Returns the name of the force field."""

        return self._forceField.name().c_str()

    ### Static Methods ########################################################
    @classmethod
    def create(cls, char *name):
        """Creates a new force field."""

        cdef _ForceField *_forceField = _ForceField_create(name)
        if _forceField is NULL:
            return None

        cdef ForceField descriptor = ForceField.__new__(ForceField)
        descriptor._forceField = _forceField

        return descriptor

    @classmethod
    def forceFields(cls):
        """Returns a list of the names of supported force fields."""

        cdef vector[string] _forceFields = _ForceField_forceFields()

        forceFields = []
        for i in range(_forceFields.size()):
            forceFields.append(_forceFields[i].c_str())

        return forceFields

