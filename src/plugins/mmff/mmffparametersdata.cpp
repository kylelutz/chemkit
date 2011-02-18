/******************************************************************************
**
** Copyright (C) 2009-2011 Kyle Lutz <kyle.r.lutz@gmail.com>
**
** This file is part of chemkit. For more information see
** <http://www.chemkit.org>.
**
** chemkit is free software: you can redistribute it and/or modify
** it under the terms of the GNU Lesser General Public License as published by
** the Free Software Foundation, either version 3 of the License, or
** (at your option) any later version.
**
** chemkit is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU Lesser General Public License for more details.
**
** You should have received a copy of the GNU Lesser General Public License
** along with chemkit. If not, see <http://www.gnu.org/licenses/>.
**
******************************************************************************/

#include "mmffparametersdata.h"

// === MmffParametersData ================================================== //
// --- Construction and Destruction ---------------------------------------- //
/// Creates a new parameters data object.
MmffParametersData::MmffParametersData()
    : vanDerWaalsParameters(MmffParameters::MaxAtomType + 1),
      partialChargeParameters(MmffParameters::MaxAtomType + 1)
{
    m_refcount.ref();
}

/// Destroys the parameters data object. This should not be called
/// directly, instead use the deref() method.
MmffParametersData::~MmffParametersData()
{
    qDeleteAll(bondStrechParameters.values());
    qDeleteAll(angleBendParameters.values());
    qDeleteAll(strechBendParameters.values());
    qDeleteAll(defaultStrechBendParameters);
    qDeleteAll(outOfPlaneBendingParameters.values());
    qDeleteAll(torsionParameters.values());
    qDeleteAll(vanDerWaalsParameters);
    qDeleteAll(chargeParameters);
    qDeleteAll(partialChargeParameters);
}

// --- Reference Counting -------------------------------------------------- //
/// Increases the reference count by one.
void MmffParametersData::ref()
{
    m_refcount.ref();
}

/// Decreases the reference count by one and deletes the object if
/// the reference count is zero.
void MmffParametersData::deref()
{
    if(!m_refcount.deref()){
        delete this;
    }
}
