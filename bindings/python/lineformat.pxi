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
from libcpp.vector cimport vector

from lineformat cimport _LineFormat
from lineformat cimport create as _LineFormat_create
from lineformat cimport formats as _LineFormat_formats

cdef class LineFormat:
    cdef _LineFormat *_lineFormat

    ### Construction and Destruction ##########################################
    def __cinit__(self):
        self._lineFormat = NULL

    def __init__(self):
        raise TypeError("LineFormat objects are not constructible, use the LineFormat.create() method")

    def __dealloc__(self):
        del self._lineFormat

    ### Properties ############################################################
    def name(self):
        """Returns the name of the line format."""

        return self._lineFormat.name().c_str()

    ### Input and Output ######################################################
    def read(self, char *formula):
        """Reads the formula and returns a new molecule."""

        return Molecule_fromPointer(self._lineFormat.read(formula))

    def write(self, Molecule molecule):
        """Returns the formula for the molecule."""

        return self._lineFormat.write(molecule._molecule).c_str()

    ### Error Handling ########################################################
    def errorString(self):
        """Returns a string describing the last error that occurred."""

        return self._lineFormat.errorString().c_str()

    ### Static Methods ########################################################
    @classmethod
    def create(cls, char *name):
        """Creates a new line format."""

        cdef _LineFormat *_lineFormat = _LineFormat_create(name)
        if _lineFormat is NULL:
            return None

        cdef LineFormat f = LineFormat.__new__(LineFormat)
        f._lineFormat = _lineFormat
        return f

    @classmethod
    def formats(cls):
        """Returns a list of the names of the supported line formats."""

        cdef vector[string] _formats = _LineFormat_formats()

        formats = []
        for i in range(_formats.size()):
            formats.append(_formats[i].c_str())

        return formats

