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

from fingerprint cimport _Fingerprint
from fingerprint cimport create as _Fingerprint_create
from fingerprint cimport fingerprints as _Fingerprint_fingerprints

cdef class Fingerprint:
    cdef _Fingerprint *_fingerprint

    ### Construction and Destruction ##########################################
    def __cinit__(self):
        self._fingerprint = NULL

    def __init__(self):
        raise TypeError("Fingerprint objects are not constructible, use the Fingerprint.create() method")

    def __dealloc__(self):
        del self._fingerprint

    ### Properties ############################################################
    def name(self):
        """Returns the name of the fingerprint."""

        return self._fingerprint.name().c_str()

    ### Static Methods ########################################################
    @classmethod
    def create(cls, char *name):
        """Creates a new fingerprint."""

        cdef _Fingerprint *_fingerprint = _Fingerprint_create(name)
        if _fingerprint is NULL:
            return None

        cdef Fingerprint descriptor = Fingerprint.__new__(Fingerprint)
        descriptor._fingerprint = _fingerprint

        return descriptor

    @classmethod
    def fingerprints(cls):
        """Returns a list of the names of supported fingerprints."""

        cdef vector[string] _fingerprints = _Fingerprint_fingerprints()

        fingerprints = []
        for i in range(_fingerprints.size()):
            fingerprints.append(_fingerprints[i].c_str())

        return fingerprints

