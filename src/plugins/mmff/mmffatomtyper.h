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

#ifndef MMFFATOMTYPER_H
#define MMFFATOMTYPER_H

#include <vector>

#include <chemkit/ring.h>
#include <chemkit/atomtyper.h>

class MmffAtomTyper : public chemkit::AtomTyper
{
public:
    // construction and destruction
    MmffAtomTyper(const chemkit::Molecule *molecule = 0);
    ~MmffAtomTyper();

    // properties
    void setMolecule(const chemkit::Molecule *molecule) CHEMKIT_OVERRIDE;

    // types
    int typeNumber(const chemkit::Atom *atom) const CHEMKIT_OVERRIDE;

    // charges
    chemkit::Real formalCharge(int index) const;
    chemkit::Real formalCharge(const chemkit::Atom *atom) const;

private:
    void setType(int index, int type, chemkit::Real formalCharge = 0);
    void setType(int index, const chemkit::Atom *atom);
    void setHydrogenType(int index, const chemkit::Atom *atom);
    void setCarbonType(int index, const chemkit::Atom *atom);
    void setNitrogenType(int index, const chemkit::Atom *atom);
    void setOxygenType(int index, const chemkit::Atom *atom);
    void setSulfurType(int index, const chemkit::Atom *atom);
    void setAromaticType(int index, const chemkit::Atom *atom, const chemkit::Ring *ring, int position);

private:
    std::vector<int> m_types;
    std::vector<chemkit::Real> m_formalCharges;
};

#endif // MMFFATOMTYPER_H
