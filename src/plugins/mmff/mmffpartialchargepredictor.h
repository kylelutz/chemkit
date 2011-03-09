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

#ifndef MMFFPARTIALCHARGEPREDICTOR_H
#define MMFFPARTIALCHARGEPREDICTOR_H

#include <chemkit/partialchargepredictor.h>

#include "mmffatomtyper.h"
#include "mmffparameters.h"

class MmffPartialChargePredictor : public chemkit::PartialChargePredictor
{
    public:
        // construction and destruction
        MmffPartialChargePredictor();
        ~MmffPartialChargePredictor();

        // properties
        void setAtomTyper(const MmffAtomTyper *typer);

        // partial charges
        virtual chemkit::Float partialCharge(int index) const;
        virtual chemkit::Float partialCharge(const chemkit::Atom *atom) const;

    protected:
        virtual void assignPartialCharges(const chemkit::Molecule *molecule);

    private:
        QVector<chemkit::Float> m_partialCharges;
        const MmffAtomTyper *m_typer;
        MmffParameters *m_parameters;
};

#endif // MMFFPARTIALCHARGEPREDICTOR_H
