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

#include "oplscalculation.h"

#include <boost/lexical_cast.hpp>

// === OplsCalculation ===================================================== //
OplsCalculation::OplsCalculation(int type, int atomCount, int parameterCount)
    : chemkit::ForceFieldCalculation(type, atomCount, parameterCount)
{
}

// === OplsBondStrechCalculation =========================================== //
OplsBondStrechCalculation::OplsBondStrechCalculation(const chemkit::ForceFieldAtom *a, const chemkit::ForceFieldAtom *b)
    : OplsCalculation(BondStrech, 2, 2)
{
    setAtom(0, a);
    setAtom(1, b);
}

bool OplsBondStrechCalculation::setup(const OplsParameters *parameters)
{
    const chemkit::ForceFieldAtom *a = atom(0);
    const chemkit::ForceFieldAtom *b = atom(1);

    int typeA = boost::lexical_cast<int>(a->type());
    int typeB = boost::lexical_cast<int>(b->type());

    const OplsBondStrechParameters *p = parameters->bondStrechParameters(typeA, typeB);
    if(!p){
        return false;
    }

    setParameter(0, p->kb);
    setParameter(1, p->r0);

    return true;
}

chemkit::Float OplsBondStrechCalculation::energy() const
{
    const chemkit::ForceFieldAtom *a = atom(0);
    const chemkit::ForceFieldAtom *b = atom(1);

    chemkit::Float kb = parameter(0);
    chemkit::Float r0 = parameter(1);

    chemkit::Float r = distance(a, b);

    return kb * pow(r - r0, 2);
}

QVector<chemkit::Vector3> OplsBondStrechCalculation::gradient() const
{
    const chemkit::ForceFieldAtom *a = atom(0);
    const chemkit::ForceFieldAtom *b = atom(1);

    chemkit::Float kb = parameter(0);
    chemkit::Float r0 = parameter(1);

    chemkit::Float r = distance(a, b);

    QVector<chemkit::Vector3> gradient = distanceGradient(a, b);

    // dE/dr
    chemkit::Float de_dr = 2.0 * kb * (r - r0);

    gradient[0] *= de_dr;
    gradient[1] *= de_dr;

    return gradient;
}

// === OplsAngleBendCalculation ============================================ //
OplsAngleBendCalculation::OplsAngleBendCalculation(const chemkit::ForceFieldAtom *a, const chemkit::ForceFieldAtom *b, const chemkit::ForceFieldAtom *c)
    : OplsCalculation(AngleBend, 3, 2)
{
    setAtom(0, a);
    setAtom(1, b);
    setAtom(2, c);
}

bool OplsAngleBendCalculation::setup(const OplsParameters *parameters)
{
    const chemkit::ForceFieldAtom *a = atom(0);
    const chemkit::ForceFieldAtom *b = atom(1);
    const chemkit::ForceFieldAtom *c = atom(2);

    int typeA = boost::lexical_cast<int>(a->type());
    int typeB = boost::lexical_cast<int>(b->type());
    int typeC = boost::lexical_cast<int>(c->type());

    const OplsAngleBendParameters *p = parameters->angleBendParameters(typeA, typeB, typeC);
    if(!p){
        return false;
    }

    setParameter(0, p->ka);
    setParameter(1, p->theta0 * chemkit::constants::DegreesToRadians);

    return true;
}

chemkit::Float OplsAngleBendCalculation::energy() const
{
    const chemkit::ForceFieldAtom *a = atom(0);
    const chemkit::ForceFieldAtom *b = atom(1);
    const chemkit::ForceFieldAtom *c = atom(2);

    chemkit::Float ka = parameter(0);
    chemkit::Float theta0 = parameter(1);

    chemkit::Float theta = bondAngleRadians(a, b, c);

    return ka * pow(theta - theta0, 2);
}

QVector<chemkit::Vector3> OplsAngleBendCalculation::gradient() const
{
    const chemkit::ForceFieldAtom *a = atom(0);
    const chemkit::ForceFieldAtom *b = atom(1);
    const chemkit::ForceFieldAtom *c = atom(2);

    chemkit::Float ka = parameter(0);
    chemkit::Float theta0 = parameter(1);

    chemkit::Float theta = bondAngleRadians(a, b, c);

    QVector<chemkit::Vector3> gradient = bondAngleGradientRadians(a, b, c);

    // dE/dtheta
    chemkit::Float de_dtheta = (2.0 * ka * (theta - theta0));

    gradient[0] *= de_dtheta;
    gradient[1] *= de_dtheta;
    gradient[2] *= de_dtheta;

    return gradient;
}

// === OplsTorsionCalculation ============================================== //
OplsTorsionCalculation::OplsTorsionCalculation(const chemkit::ForceFieldAtom *a, const chemkit::ForceFieldAtom *b, const chemkit::ForceFieldAtom *c, const chemkit::ForceFieldAtom *d)
    : OplsCalculation(Torsion, 4, 3)
{
    setAtom(0, a);
    setAtom(1, b);
    setAtom(2, c);
    setAtom(3, d);
}

bool OplsTorsionCalculation::setup(const OplsParameters *parameters)
{
    const chemkit::ForceFieldAtom *a = atom(0);
    const chemkit::ForceFieldAtom *b = atom(1);
    const chemkit::ForceFieldAtom *c = atom(2);
    const chemkit::ForceFieldAtom *d = atom(3);

    int typeA = boost::lexical_cast<int>(a->type());
    int typeB = boost::lexical_cast<int>(b->type());
    int typeC = boost::lexical_cast<int>(c->type());
    int typeD = boost::lexical_cast<int>(d->type());

    const OplsTorsionParameters *p = parameters->torsionParameters(typeA, typeB, typeC, typeD);
    if(!p){
        return false;
    }

    setParameter(0, p->v1);
    setParameter(1, p->v2);
    setParameter(2, p->v3);

    return true;
}

chemkit::Float OplsTorsionCalculation::energy() const
{
    const chemkit::ForceFieldAtom *a = atom(0);
    const chemkit::ForceFieldAtom *b = atom(1);
    const chemkit::ForceFieldAtom *c = atom(2);
    const chemkit::ForceFieldAtom *d = atom(3);

    chemkit::Float v1 = parameter(0);
    chemkit::Float v2 = parameter(1);
    chemkit::Float v3 = parameter(2);

    chemkit::Float phi = torsionAngleRadians(a, b, c, d);

    return (1.0/2.0) * (v1 * (1.0 + cos(phi)) + v2 * (1.0 - cos(2.0 * phi)) + v3 * (1.0 + cos(3.0 * phi)));
}

QVector<chemkit::Vector3> OplsTorsionCalculation::gradient() const
{
    const chemkit::ForceFieldAtom *a = atom(0);
    const chemkit::ForceFieldAtom *b = atom(1);
    const chemkit::ForceFieldAtom *c = atom(2);
    const chemkit::ForceFieldAtom *d = atom(3);

    chemkit::Float v1 = parameter(0);
    chemkit::Float v2 = parameter(1);
    chemkit::Float v3 = parameter(2);

    chemkit::Float phi = torsionAngleRadians(a, b, c, d);

    // dE/dphi
    chemkit::Float de_dphi = (1.0/2.0) * (-v1 * sin(phi) + 2.0 * v2 * sin(2.0 * phi) - 3.0 * v3 * sin(3.0 * phi));

    QVector<chemkit::Vector3> gradient = torsionAngleGradientRadians(a, b, c, d);

    gradient[0] *= de_dphi;
    gradient[1] *= de_dphi;
    gradient[2] *= de_dphi;
    gradient[3] *= de_dphi;

    return gradient;
}

// === OplsNonbondedCalculation ============================================ //
OplsNonbondedCalculation::OplsNonbondedCalculation(const chemkit::ForceFieldAtom *a, const chemkit::ForceFieldAtom *b)
    : OplsCalculation(VanDerWaals | Electrostatic, 2, 5)
{
    setAtom(0, a);
    setAtom(1, b);
}

bool OplsNonbondedCalculation::setup(const OplsParameters *parameters)
{
    const chemkit::ForceFieldAtom *a = atom(0);
    const chemkit::ForceFieldAtom *b = atom(1);

    int typeA = boost::lexical_cast<int>(a->type());
    int typeB = boost::lexical_cast<int>(b->type());

    const OplsVanDerWaalsParameters *pa = parameters->vanDerWaalsParameters(typeA);
    const OplsVanDerWaalsParameters *pb = parameters->vanDerWaalsParameters(typeB);
    if(!pa || !pb){
        return false;
    }

    chemkit::Float qa = parameters->partialCharge(typeA);
    chemkit::Float qb = parameters->partialCharge(typeB);
    chemkit::Float sigma = sqrt(pa->sigma * pb->sigma);
    chemkit::Float epsilon = sqrt(pa->epsilon * pb->epsilon);

    setParameter(0, qa);
    setParameter(1, qb);
    setParameter(2, sigma);
    setParameter(3, epsilon);

    // one-four scaling
    if(a->isOneFour(b)){
        setParameter(4, 0.5);
    }
    else{
        setParameter(4, 1.0);
    }

    return true;
}

chemkit::Float OplsNonbondedCalculation::energy() const
{
    const chemkit::ForceFieldAtom *a = atom(0);
    const chemkit::ForceFieldAtom *b = atom(1);

    chemkit::Float qa = parameter(0);
    chemkit::Float qb = parameter(1);
    chemkit::Float e = 332.06; // vacuum permitivity
    chemkit::Float sigma = parameter(2);
    chemkit::Float epsilon = parameter(3);
    chemkit::Float scale = parameter(4);

    chemkit::Float r = distance(a, b);

    return scale * ((qa * qb * e) / r + 4.0 * epsilon * (pow(sigma / r, 12) - pow(sigma / r, 6)));
}

QVector<chemkit::Vector3> OplsNonbondedCalculation::gradient() const
{
    QVector<chemkit::Vector3> gradient(2);

    const chemkit::ForceFieldAtom *a = atom(0);
    const chemkit::ForceFieldAtom *b = atom(1);

    chemkit::Float qa = parameter(0);
    chemkit::Float qb = parameter(1);
    chemkit::Float e = 332.06; // vacuum permitivity
    chemkit::Float sigma = parameter(2);
    chemkit::Float epsilon = parameter(3);
    chemkit::Float scale = parameter(4);

    chemkit::Float r = distance(a, b);
    chemkit::Float sr = sigma / r;

    // dE/dr
    chemkit::Float de_dr = scale * ((1.0 / pow(r, 3)) * (-qa * qb * e + -4.0 * epsilon * sigma * (12.0 * pow(sr, 11) - 6.0 * pow(sr, 5))));

    // dE/da
    chemkit::Vector3 de_da = (a->position() - b->position()) * de_dr;

    gradient[0] = de_da;
    gradient[1] = -de_da;

    return gradient;
}
