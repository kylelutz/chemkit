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

from moleculegeometryoptimizer cimport _MoleculeGeometryOptimizer
from moleculegeometryoptimizer cimport optimizeCoordinates as _MoleculeGeometryOptimizer_optimizeCoordinates

cdef class MoleculeGeometryOptimizer:
    cdef _MoleculeGeometryOptimizer *_optimizer

    ### Construction and Destruction ##########################################
    def __cinit__(self):
        self._optimizer = NULL

    def __init__(self):
        """Creates a new geometry optimizer."""

        self._optimizer = new _MoleculeGeometryOptimizer()

    def __dealloc__(self):
        """Destroys the geometry optimizer."""

        del self._optimizer

    ### Properties ############################################################
    def setMolecule(self, Molecule molecule):
        """Sets the molecule to be optimized."""

        self._optimizer.setMolecule(molecule._molecule)

    def molecule(self):
        """Returns the molecule."""

        return Molecule_fromPointer(self._optimizer.molecule())

    def setForceField(self, char *forceField):
        """Sets the force field to use for the optimization."""

        return self._optimizer.setForceField(forceField)

    def forceField(self):
        """Returns the name of the force field."""

        cdef string name = self._optimizer.forceField()

        return name.c_str()

    ### Energy ################################################################
    def energy(self):
        """Returns the current energy of the molecule."""

        return self._optimizer.energy()

    ### Optimization ##########################################################
    def setup(self):
        """Sets up the force field for the optimization."""

        return self._optimizer.setup()

    def step(self):
        """Performs a single step of optimization."""

        self._optimizer.step()

    def converged(self):
        """Returns True if the optimization algorithm has converged."""

        return self._optimizer.converged()

    def optimize(self):
        """Runs the optimization algorithm until converged."""

        return self._optimizer.optimize()

    def writeCoordinates(self):
        """Writes the optimized coordinates to the molecule."""

        self._optimizer.writeCoordinates()


    ### Error Handling ########################################################
    def errorString(self):
        """Returns a string describing the last error that occurred."""

        cdef string error = self._optimizer.errorString()

        return error.c_str()

    ### Static Methods ########################################################
    @classmethod
    def optimizeCoordinates(cls, Molecule molecule):
        """Optimizes the 3D coordinates for the atoms in molecule."""

        _MoleculeGeometryOptimizer_optimizeCoordinates(molecule._molecule)
