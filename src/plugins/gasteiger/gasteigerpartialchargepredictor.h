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

#include <QtCore>

#include <chemkit/partialchargepredictor.h>

struct GasteigerParameters {
    chemkit::Float a;
    chemkit::Float b;
    chemkit::Float c;
};

class GasteigerPartialChargePredictor : public chemkit::PartialChargePredictor
{
    public:
        // construction and destruction
        GasteigerPartialChargePredictor();
        ~GasteigerPartialChargePredictor();

        // partial charges
        chemkit::Float partialCharge(int index) const;

    protected:
        void assignPartialCharges(const chemkit::Molecule *molecule);

    private:
        const GasteigerParameters* atomParameters(const chemkit::Atom *atom) const;

    private:
        QVector<chemkit::Float> m_charges;
        QVector<chemkit::Float> m_electronegativies;
        QVector<const GasteigerParameters *> m_parameters;
};
