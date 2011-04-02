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

#ifndef CHEMKIT_BONDPREDICTOR_H
#define CHEMKIT_BONDPREDICTOR_H

#include "chemkit.h"

#include <vector>

namespace chemkit {

class Atom;
class Molecule;
class BondPredictorPrivate;

class CHEMKIT_EXPORT BondPredictor
{
    public:
        // construction and destruction
        BondPredictor(Molecule *molecule);
        ~BondPredictor();

        // properties
        void setTolerance(Float tolerance);
        Float tolerance() const;
        void setMinimumBondLength(Float length);
        Float minimumBondLength() const;
        void setMaximumBondLength(Float length);
        Float maximumBondLength() const;
        Molecule* molecule() const;

        // prediction
        std::vector<std::pair<Atom *, Atom *> > predictedBonds();
        bool couldBeBonded(Atom *a, Atom *b) const;

        // static methods
        static void predictBonds(Molecule *molecule);

    private:
        BondPredictorPrivate* const d;
};

} // end chemkit namespace

#endif // CHEMKIT_BONDPREDICTOR_H
