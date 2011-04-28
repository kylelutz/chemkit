/******************************************************************************
**
** Copyright (C) 2009-2011 Kyle Lutz <kyle.r.lutz@gmail.com>
** All rights reserved.
**
** This file is a part of the chemkit project. For more information
** see <http://www.chemkit.org>.
**
** Redistribution and use in source and binary forms, with or without
** modification, are permitted provided that the following conditions
** are met:
**
**   * Redistributions of source code must retain the above copyright
**     notice, this list of conditions and the following disclaimer.
**   * Redistributions in binary form must reproduce the above copyright
**     notice, this list of conditions and the following disclaimer in the
**     documentation and/or other materials provided with the distribution.
**   * Neither the name of the chemkit project nor the names of its
**     contributors may be used to endorse or promote products derived
**     from this software without specific prior written permission.
**
** THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
** "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
** LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
** A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
** OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
** SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
** LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
** DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
** THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
** (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
** OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
**
******************************************************************************/

#include "unitcell.h"

namespace chemkit {

// === UnitCellPrivate ===================================================== //
class UnitCellPrivate
{
    public:
        Vector3 x;
        Vector3 y;
        Vector3 z;
};

// === UnitCell ============================================================ //
/// \class UnitCell unitcell.h chemkit/unitcell.h
/// \ingroup chemkit
/// \brief The UnitCell class represents a unit cell.

// --- Construction and Destruction ---------------------------------------- //
/// Creates a new unit cell.
UnitCell::UnitCell()
    : d(new UnitCellPrivate)
{
}

/// Creates a new unit cell with \p x, \p y, and \p z.
UnitCell::UnitCell(const Vector3 &x, const Vector3 &y, const Vector3 &z)
    : d(new UnitCellPrivate)
{
    d->x = x;
    d->y = y;
    d->z = z;
}

/// Destroys the unit cell object.
UnitCell::~UnitCell()
{
    delete d;
}

// --- Properties ---------------------------------------------------------- //
/// Returns the x-vector for the unit cell.
const Vector3& UnitCell::x() const
{
    return d->x;
}

/// Returns the y-vector for the unit cell.
const Vector3& UnitCell::y() const
{
    return d->y;
}

/// Returns the z-vector for the unit cell.
const Vector3& UnitCell::z() const
{
    return d->z;
}

} // end chemkit namespace
