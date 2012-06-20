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

from partialchargemodel cimport _PartialChargeModel
from partialchargemodel cimport create as _PartialChargeModel_create
from partialchargemodel cimport models as _PartialChargeModel_models
from partialchargemodel cimport assignPartialCharges as _PartialChargeModel_assignPartialCharges

cdef class PartialChargeModel:
    """The PartialChargeModel class assigns partial charges to atoms."""

    cdef _PartialChargeModel *_model

    ### Construction and Destruction ##########################################
    def __cinit__(self):
        self._model = NULL

    def __init__(self):
        raise TypeError("PartialChargeModel objects are not constructible, "
                        "use the PartialChargeModel.create() method")

    def __dealloc__(self):
        del self._model

    ### Properties ############################################################
    def name(self):
        """Returns the name of the partial charge model."""

        cdef string name = self._model.name()

        return name.c_str()

    def setMolecule(self, Molecule molecule):
        """Sets the molecule."""

        self._model.setMolecule(molecule._molecule)

    def molecule(self):
        """Returns the molecule."""

        return Molecule_fromPointer(self._model.molecule())

    ### Partial Charges #######################################################
    def partialCharge(self, Atom atom):
        """Returns the partial charge for the atom."""

        return self._model.partialCharge(atom._atom)

    ### Static Methods ########################################################
    @classmethod
    def create(cls, char *name):
        """Creates a new atom typer."""

        cdef _PartialChargeModel *_model = _PartialChargeModel_create(name)
        if _model is NULL:
            return None

        cdef PartialChargeModel model = PartialChargeModel.__new__(PartialChargeModel)
        model._model = _model
        return model

    @classmethod
    def models(cls):
        """Returns a list of the names of the supported partial charge models."""

        cdef vector[string] _models = _PartialChargeModel_models()

        models = []
        for i in range(_models.size()):
            models.append(_models[i].c_str())

        return models

    @classmethod
    def assignPartialCharges(cls, Molecule molecule, char *model):
        """Assigns partial charges to the atoms in the molecule using model."""

        return _PartialChargeModel_assignPartialCharges(molecule._molecule, <string>(model))
