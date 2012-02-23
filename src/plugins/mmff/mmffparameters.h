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

#ifndef MMFFPARAMETERS_H
#define MMFFPARAMETERS_H

#include <string>

#include "mmffatom.h"
#include "mmffaromaticitymodel.h"

class MmffParametersData;

struct MmffBondStrechParameters
{
    chemkit::Real kb;
    chemkit::Real r0;
};

struct MmffAngleBendParameters
{
    chemkit::Real ka;
    chemkit::Real theta0;
};

struct MmffStrechBendParameters
{
    chemkit::Real kba_ijk;
    chemkit::Real kba_kji;
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
    chemkit::Real koop;
};

struct MmffTorsionParameters
{
    chemkit::Real V1;
    chemkit::Real V2;
    chemkit::Real V3;
};

struct MmffVanDerWaalsParameters
{
    chemkit::Real alpha;
    chemkit::Real N;
    chemkit::Real A;
    chemkit::Real G;
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
    chemkit::Real bci;
};

struct MmffPartialChargeParameters
{
    chemkit::Real pbci;
    chemkit::Real fcadj;
};

class MmffParameters
{
public:
    // construction and destruction
    MmffParameters();
    ~MmffParameters();

    // parameters
    std::string fileName() const;
    bool read(const std::string &fileName);
    const MmffBondStrechParameters* bondStrechParameters(const MmffAtom *a, const MmffAtom *b) const;
    const MmffAngleBendParameters* angleBendParameters(const MmffAtom *a, const MmffAtom *b, const MmffAtom *c) const;
    const MmffStrechBendParameters* strechBendParameters(const MmffAtom *a, const MmffAtom *b, const MmffAtom *c) const;
    const MmffStrechBendParameters* defaultStrechBendParameters(const MmffAtom *a, const MmffAtom *b, const MmffAtom *c) const;
    const MmffOutOfPlaneBendingParameters* outOfPlaneBendingParameters(const MmffAtom *a, const MmffAtom *b, const MmffAtom *c, const MmffAtom *d) const;
    const MmffTorsionParameters* torsionParameters(const MmffAtom *a, const MmffAtom *b, const MmffAtom *c, const MmffAtom *d) const;
    const MmffVanDerWaalsParameters* vanDerWaalsParameters(const MmffAtom *atom) const;
    const MmffAtomParameters* atomParameters(int type) const;
    const MmffAtomParameters* atomParameters(const MmffAtom *atom) const;
    const MmffChargeParameters* chargeParameters(const chemkit::Atom *a, int typeA, const chemkit::Atom *b, int typeB) const;
    const MmffChargeParameters* chargeParameters(const MmffAtom *a, const MmffAtom *b) const;
    const MmffPartialChargeParameters* partialChargeParameters(int type) const;
    const MmffPartialChargeParameters* partialChargeParameters(const MmffAtom *atom) const;

    // error handling
    std::string errorString() const;

    // constants
    const static int MaxAtomType = 99;

private:
    const MmffBondStrechParameters* bondStrechParameters(int bondType, int typeA, int typeB) const;
    const MmffBondStrechParameters* empiricalBondStrechParameters(int atomicNumberA, int atomicNumberB) const;
    const MmffAngleBendParameters* angleBendParameters(int angleType, int typeA, int typeB, int typeC) const;
    const MmffStrechBendParameters* strechBendParameters(int strechBendType, int typeA, int typeB, int typeC) const;
    const MmffStrechBendParameters* defaultStrechBendParameters(int rowA, int rowB, int rowC) const;
    const MmffOutOfPlaneBendingParameters* outOfPlaneBendingParameters(int typeA, int typeB, int typeC, int typeD) const;
    const MmffTorsionParameters* torsionParameters(int torsionType, int typeA, int typeB, int typeC, int typeD) const;
    int calculateBondType(const chemkit::Bond *bond, int typeA, int typeB) const;
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
    void setErrorString(const std::string &errorString);

private:
    std::string m_fileName;
    std::string m_errorString;
    boost::shared_ptr<MmffParametersData> d;
    MmffAromaticityModel m_aromaticityModel;
};

#endif // MMFFPARAMETERS_H
