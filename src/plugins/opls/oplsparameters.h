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

#ifndef OPLSPARAMETERS_H
#define OPLSPARAMETERS_H

#include <QtCore>

#include <chemkit/chemkit.h>
#include <chemkit/forcefieldatom.h>

struct OplsBondStrechParameters
{
    int typeA;
    int typeB;
    chemkit::Float kb;
    chemkit::Float r0;
};

struct OplsAngleBendParameters
{
    int typeA;
    int typeB;
    int typeC;
    chemkit::Float ka;
    chemkit::Float theta0;
};

struct OplsTorsionParameters
{
    int typeA;
    int typeB;
    int typeC;
    int typeD;
    chemkit::Float v1;
    chemkit::Float v2;
    chemkit::Float v3;
};

struct OplsVanDerWaalsParameters
{
    chemkit::Float sigma;
    chemkit::Float epsilon;
};

class OplsParameters
{
    public:
        // construction and destruction
        OplsParameters(const QString &fileName);
        ~OplsParameters();

        // properties
        void setFileName(const QString &fileName);
        QString fileName() const;

        // parameters
        int atomClass(int type) const;
        QString atomName(int type) const;
        chemkit::Float partialCharge(int type) const;
        const OplsBondStrechParameters* bondStrechParameters(int a, int b) const;
        const OplsAngleBendParameters* angleBendParameters(int a, int b, int c) const;
        const OplsTorsionParameters* torsionParameters(int a, int b, int c, int d) const;
        const OplsVanDerWaalsParameters* vanDerWaalsParameters(int type) const;

    private:
        bool read(const QString &fileName);

    private:
        QString m_fileName;
        QVector<int> m_typeToClass;
        QVector<QString> m_typeToName;
        QVector<OplsBondStrechParameters> m_bondStrechParameters;
        QVector<OplsAngleBendParameters> m_angleBendParameters;
        QVector<OplsTorsionParameters> m_torsionParameters;
        QVector<OplsVanDerWaalsParameters> m_vanDerWaalsParameters;
        QVector<chemkit::Float> m_typeToCharge;
};

#endif // OPLSPARAMETERS_H
