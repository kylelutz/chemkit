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

#include "mmffatom.h"

#include <boost/lexical_cast.hpp>

#include <chemkit/atom.h>

#include "mmffforcefield.h"
#include "mmffparameters.h"

// --- Construction and Destruction ---------------------------------------- //
MmffAtom::MmffAtom(chemkit::md::ForceField *forceField, const chemkit::Atom *atom)
    : chemkit::md::ForceFieldAtom(forceField, atom)
{
    m_typeNumber = 0;
    m_formalCharge = 0;
}

// --- Properties ---------------------------------------------------------- //
const MmffForceField* MmffAtom::forceField() const
{
    return static_cast<const MmffForceField *>(chemkit::md::ForceFieldAtom::forceField());
}

void MmffAtom::setType(int typeNumber, chemkit::Real formalCharge)
{
    m_typeNumber = typeNumber;
    m_formalCharge = formalCharge;
}

std::string MmffAtom::type() const
{
    return boost::lexical_cast<std::string>(m_typeNumber);
}

int MmffAtom::typeNumber() const
{
    return m_typeNumber;
}

chemkit::Real MmffAtom::formalCharge() const
{
    return m_formalCharge;
}

int MmffAtom::period() const
{
    return atom()->element().period();
}

// --- Parameters ---------------------------------------------------------- //
const MmffAtomParameters* MmffAtom::parameters() const
{
    const MmffParameters *mmffParameters = static_cast<const MmffForceField *>(forceField())->parameters();
    if(!mmffParameters)
        return 0;

    return mmffParameters->atomParameters(this);
}
