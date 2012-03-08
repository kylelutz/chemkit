/******************************************************************************
**
** Copyright (C) 2009-2012 Kyle Lutz <kyle.r.lutz@gmail.com>
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

#include "vabcdescriptor.h"

#include <boost/bind.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/range/algorithm.hpp>

#include <chemkit/atom.h>
#include <chemkit/bond.h>
#include <chemkit/ring.h>
#include <chemkit/foreach.h>
#include <chemkit/molecule.h>
#include <chemkit/aromaticitymodel.h>

VabcDescriptor::VabcDescriptor()
    : chemkit::MolecularDescriptor("vabc")
{
    // values from table 2 and the supplementary spreadsheet
    m_volumes[chemkit::Atom::Hydrogen] = 7.2382293504;
    m_volumes[chemkit::Atom::Carbon] = 20.5795259250667;
    m_volumes[chemkit::Atom::Nitrogen] = 15.5985308577667;
    m_volumes[chemkit::Atom::Oxygen] = 14.7102267005611;
    m_volumes[chemkit::Atom::Fluorine] = 13.3057882007064;
    m_volumes[chemkit::Atom::Chlorine] = 22.4492971208333;
    m_volumes[chemkit::Atom::Bromine] = 26.5218483279667;
    m_volumes[chemkit::Atom::Iodine] = 32.5150310206656;
    m_volumes[chemkit::Atom::Phosphorus] = 24.4290240576;
    m_volumes[chemkit::Atom::Sulfur] = 24.4290240576;
    m_volumes[chemkit::Atom::Arsenic] = 26.5218483279667;
    m_volumes[chemkit::Atom::Boron] = 40.48;
    m_volumes[chemkit::Atom::Silicon] = 38.7923854248;
    m_volumes[chemkit::Atom::Selenium] = 28.7309115245333;
    m_volumes[chemkit::Atom::Tellurium] = 36.62;
}

VabcDescriptor::~VabcDescriptor()
{
}

chemkit::Real VabcDescriptor::volume(const chemkit::Atom *atom) const
{
    std::map<chemkit::Element::AtomicNumberType, chemkit::Real>::const_iterator iter =
        m_volumes.find(atom->atomicNumber());

    if(iter != m_volumes.end()){
        return iter->second;
    }
    else{
        return 0;
    }
}

// Returns the VABC value for the molecule.
chemkit::Variant VabcDescriptor::value(const chemkit::Molecule *molecule) const
{
    chemkit::Real atomContributions = 0;

    foreach(const chemkit::Atom *atom, molecule->atoms()){
        atomContributions += volume(atom);
    }

    size_t bondCount = molecule->bondCount();
    size_t ringCount = molecule->ringCount();
    size_t aromaticRingCount = boost::count_if(molecule->rings(), boost::bind(&chemkit::Ring::isAromatic, _1));
    size_t aliphaticRingCount = ringCount - aromaticRingCount;

    // equation 4
    return atomContributions - (5.92 * bondCount) - (14.7 * aromaticRingCount) - (3.8 * aliphaticRingCount);
}
