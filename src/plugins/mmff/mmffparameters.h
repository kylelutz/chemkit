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

#ifndef MMFFPARAMETERS_H
#define MMFFPARAMETERS_H

#include <QtCore>

#include "mmffatom.h"

struct MmffBondStrechParameters
{
    chemkit::Float kb;
    chemkit::Float r0;
};

struct MmffAngleBendParameters
{
    chemkit::Float ka;
    chemkit::Float theta0;
};

struct MmffStrechBendParameters
{
    chemkit::Float kba_ijk;
    chemkit::Float kba_kji;
};

struct MmffDefaultStrechBendParameters
{
    int rowA;
    int rowB;
    int rowC;
    MmffStrechBendParameters parameters;
};

struct MmffOutOfPlaneBendingParameters
{
    chemkit::Float koop;
};

struct MmffTorsionParameters
{
    chemkit::Float V1;
    chemkit::Float V2;
    chemkit::Float V3;
};

struct MmffVanDerWaalsParameters
{
    chemkit::Float alpha;
    chemkit::Float N;
    chemkit::Float A;
    chemkit::Float G;
    char DA;
};

struct MmffAtomParameters
{
    int aspec;
    int crd;
    int val;
    int pilp;
    int mltb;
    int arom;
    int lin;
    int sbmb;
};

struct MmffChargeParameters
{
    int bondType;
    int typeA;
    int typeB;
    chemkit::Float bci;
};

struct MmffPartialChargeParameters
{
    chemkit::Float pbci;
    chemkit::Float fcadj;
};

class MmffParameters
{
    public:
        // construction and destruction
        MmffParameters();
        ~MmffParameters();

        // parameters
        QString fileName() const;
        bool read(const QString &fileName);
        const MmffBondStrechParameters* bondStrechParameters(const MmffAtom *a, const MmffAtom *b) const;
        const MmffAngleBendParameters* angleBendParameters(const MmffAtom *a, const MmffAtom *b, const MmffAtom *c) const;
        const MmffStrechBendParameters* strechBendParameters(const MmffAtom *a, const MmffAtom *b, const MmffAtom *c) const;
        const MmffStrechBendParameters* defaultStrechBendParameters(const MmffAtom *a, const MmffAtom *b, const MmffAtom *c) const;
        const MmffOutOfPlaneBendingParameters* outOfPlaneBendingParameters(const MmffAtom *a, const MmffAtom *b, const MmffAtom *c, const MmffAtom *d) const;
        const MmffTorsionParameters* torsionParameters(const MmffAtom *a, const MmffAtom *b, const MmffAtom *c, const MmffAtom *d) const;
        const MmffVanDerWaalsParameters* vanDerWaalsParameters(const MmffAtom *atom) const;
        const MmffAtomParameters* atomParameters(const MmffAtom *atom) const;
        const MmffChargeParameters* chargeParameters(const MmffAtom *a, const MmffAtom *b) const;
        const MmffPartialChargeParameters* partialChargeParameters(const MmffAtom *atom) const;

        // error handling
        QString errorString() const;

    private:
        const MmffBondStrechParameters* bondStrechParameters(int bondType, int typeA, int typeB) const;
        const MmffBondStrechParameters* empiricalBondStrechParameters(int atomicNumberA, int atomicNumberB) const;
        const MmffAngleBendParameters* angleBendParameters(int angleType, int typeA, int typeB, int typeC) const;
        const MmffStrechBendParameters* strechBendParameters(int strechBendType, int typeA, int typeB, int typeC) const;
        const MmffStrechBendParameters* defaultStrechBendParameters(int rowA, int rowB, int rowC) const;
        const MmffOutOfPlaneBendingParameters* outOfPlaneBendingParameters(int typeA, int typeB, int typeC, int typeD) const;
        const MmffTorsionParameters* torsionParameters(int torsionType, int typeA, int typeB, int typeC, int typeD) const;
        int calculateBondType(const MmffAtom *a, const MmffAtom *b) const;
        int calculateAngleType(const MmffAtom *a, const MmffAtom *b, const MmffAtom *c) const;
        int calculateStrechBendType(const MmffAtom *a, const MmffAtom *b, const MmffAtom *c) const;
        int calculateTorsionType(const MmffAtom *a, const MmffAtom *b, const MmffAtom *c, const MmffAtom *d) const;
        int equivalentType(const MmffAtom *atom, int level) const;
        int calculateBondStrechIndex(int bondType, int typeA, int typeB) const;
        int calculateAngleBendIndex(int angleType, int typeA, int typeB, int typeC) const;
        int calculateStrechBendIndex(int strechBendType, int typeA, int typeB, int typeC) const;
        int calculateOutOfPlaneBendingIndex(int typeA, int typeB, int typeC, int typeD) const;
        int calculateTorsionIndex(int torsionType, int typeA, int typeB, int typeC, int typeD) const;
        void setErrorString(const QString &errorString);

    private:
        QString m_fileName;
        QString m_errorString;
        QMap<int, MmffBondStrechParameters *> m_bondStrechParameters;
        QMap<int, MmffAngleBendParameters *> m_angleBendParameters;
        QMap<int, MmffStrechBendParameters *> m_strechBendParameters;
        QList<MmffDefaultStrechBendParameters *> m_defaultStrechBendParameters;
        QMap<int, MmffOutOfPlaneBendingParameters *> m_outOfPlaneBendingParameters;
        QMap<int, MmffTorsionParameters *> m_torsionParameters;
        QVector<MmffVanDerWaalsParameters *> m_vanDerWaalsParameters;
        QList<MmffChargeParameters *> m_chargeParameters;
        QVector<MmffPartialChargeParameters *> m_partialChargeParameters;

        const static int MaxAtomType = 99;
};

#endif // MMFFPARAMETERS_H
