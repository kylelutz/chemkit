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

#ifndef CHEMKIT_PARTIALCHARGEPREDICTOR_H
#define CHEMKIT_PARTIALCHARGEPREDICTOR_H

#include "chemkit.h"

#include <string>
#include <vector>

namespace chemkit {

class Atom;
class Molecule;
class PartialChargePredictorPrivate;

class CHEMKIT_EXPORT PartialChargePredictor
{
    public:
        // typedefs
        typedef PartialChargePredictor* (*CreateFunction)();

        // construction and destruction
        virtual ~PartialChargePredictor();

        // properties
        std::string name() const;
        void setMolecule(const Molecule *molecule);
        const Molecule* molecule() const;

        // partial charges
        virtual Float partialCharge(int index) const;
        virtual Float partialCharge(const Atom *atom) const;

        // static methods
        static PartialChargePredictor* create(const std::string &name);
        static std::vector<std::string> predictors();
        static bool predictPartialCharges(Molecule *molecule, const std::string &predictorName);

    protected:
        PartialChargePredictor(const std::string &name);
        virtual void assignPartialCharges(const Molecule *molecule);

    private:
        PartialChargePredictorPrivate* const d;
};

} // end chemkit namespace

#endif // CHEMKIT_PARTIALCHARGEPREDICTOR_H
