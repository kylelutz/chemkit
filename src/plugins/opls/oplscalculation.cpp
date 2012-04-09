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

#include "oplscalculation.h"

#include <boost/lexical_cast.hpp>

#include <chemkit/topology.h>
#include <chemkit/constants.h>
#include <chemkit/cartesiancoordinates.h>

// === OplsCalculation ===================================================== //
OplsCalculation::OplsCalculation(int type, int atomCount, int parameterCount)
    : chemkit::ForceFieldCalculation(type, atomCount, parameterCount)
{
}

// === OplsBondStrechCalculation =========================================== //
OplsBondStrechCalculation::OplsBondStrechCalculation(size_t a, size_t b)
    : OplsCalculation(BondStrech, 2, 2)
{
    setAtom(0, a);
    setAtom(1, b);
}

bool OplsBondStrechCalculation::setup(const OplsParameters *parameters)
{
    int typeA = boost::lexical_cast<int>(atomType(0));
    int typeB = boost::lexical_cast<int>(atomType(1));

    const OplsBondStrechParameters *p = parameters->bondStrechParameters(typeA, typeB);
    if(!p){
        return false;
    }

    setParameter(0, p->kb);
    setParameter(1, p->r0);

    return true;
}

chemkit::Real OplsBondStrechCalculation::energy(const chemkit::CartesianCoordinates *coordinates) const
{
    size_t a = atom(0);
    size_t b = atom(1);

    chemkit::Real kb = parameter(0);
    chemkit::Real r0 = parameter(1);

    chemkit::Real r = coordinates->distance(a, b);

    return kb * pow(r - r0, 2);
}

std::vector<chemkit::Vector3> OplsBondStrechCalculation::gradient(const chemkit::CartesianCoordinates *coordinates) const
{
    size_t a = atom(0);
    size_t b = atom(1);

    chemkit::Real kb = parameter(0);
    chemkit::Real r0 = parameter(1);

    chemkit::Real r = coordinates->distance(a, b);

    boost::array<chemkit::Vector3, 2> gradient = coordinates->distanceGradient(a, b);

    // dE/dr
    chemkit::Real de_dr = 2.0 * kb * (r - r0);

    gradient[0] *= de_dr;
    gradient[1] *= de_dr;

    return std::vector<chemkit::Vector3>(gradient.begin(), gradient.end());
}

// === OplsAngleBendCalculation ============================================ //
OplsAngleBendCalculation::OplsAngleBendCalculation(size_t a, size_t b, size_t c)
    : OplsCalculation(AngleBend, 3, 2)
{
    setAtom(0, a);
    setAtom(1, b);
    setAtom(2, c);
}

bool OplsAngleBendCalculation::setup(const OplsParameters *parameters)
{
    int typeA = boost::lexical_cast<int>(atomType(0));
    int typeB = boost::lexical_cast<int>(atomType(1));
    int typeC = boost::lexical_cast<int>(atomType(2));

    const OplsAngleBendParameters *p = parameters->angleBendParameters(typeA, typeB, typeC);
    if(!p){
        return false;
    }

    setParameter(0, p->ka);
    setParameter(1, p->theta0 * chemkit::constants::DegreesToRadians);

    return true;
}

chemkit::Real OplsAngleBendCalculation::energy(const chemkit::CartesianCoordinates *coordinates) const
{
    size_t a = atom(0);
    size_t b = atom(1);
    size_t c = atom(2);

    chemkit::Real ka = parameter(0);
    chemkit::Real theta0 = parameter(1);

    chemkit::Real theta = coordinates->angleRadians(a, b, c);

    return ka * pow(theta - theta0, 2);
}

std::vector<chemkit::Vector3> OplsAngleBendCalculation::gradient(const chemkit::CartesianCoordinates *coordinates) const
{
    size_t a = atom(0);
    size_t b = atom(1);
    size_t c = atom(2);

    chemkit::Real ka = parameter(0);
    chemkit::Real theta0 = parameter(1);

    chemkit::Real theta = coordinates->angleRadians(a, b, c);

    boost::array<chemkit::Vector3, 3> gradient = coordinates->angleGradientRadians(a, b, c);

    // dE/dtheta
    chemkit::Real de_dtheta = (2.0 * ka * (theta - theta0));

    gradient[0] *= de_dtheta;
    gradient[1] *= de_dtheta;
    gradient[2] *= de_dtheta;

    return std::vector<chemkit::Vector3>(gradient.begin(), gradient.end());
}

// === OplsTorsionCalculation ============================================== //
OplsTorsionCalculation::OplsTorsionCalculation(size_t a, size_t b, size_t c, size_t d)
    : OplsCalculation(Torsion, 4, 3)
{
    setAtom(0, a);
    setAtom(1, b);
    setAtom(2, c);
    setAtom(3, d);
}

bool OplsTorsionCalculation::setup(const OplsParameters *parameters)
{
    int typeA = boost::lexical_cast<int>(atomType(0));
    int typeB = boost::lexical_cast<int>(atomType(1));
    int typeC = boost::lexical_cast<int>(atomType(2));
    int typeD = boost::lexical_cast<int>(atomType(3));

    const OplsTorsionParameters *p = parameters->torsionParameters(typeA, typeB, typeC, typeD);
    if(!p){
        return false;
    }

    setParameter(0, p->v1);
    setParameter(1, p->v2);
    setParameter(2, p->v3);

    return true;
}

chemkit::Real OplsTorsionCalculation::energy(const chemkit::CartesianCoordinates *coordinates) const
{
    size_t a = atom(0);
    size_t b = atom(1);
    size_t c = atom(2);
    size_t d = atom(3);

    chemkit::Real v1 = parameter(0);
    chemkit::Real v2 = parameter(1);
    chemkit::Real v3 = parameter(2);

    chemkit::Real phi = coordinates->torsionAngleRadians(a, b, c, d);

    return (1.0/2.0) * (v1 * (1.0 + cos(phi)) + v2 * (1.0 - cos(2.0 * phi)) + v3 * (1.0 + cos(3.0 * phi)));
}

std::vector<chemkit::Vector3> OplsTorsionCalculation::gradient(const chemkit::CartesianCoordinates *coordinates) const
{
    size_t a = atom(0);
    size_t b = atom(1);
    size_t c = atom(2);
    size_t d = atom(3);

    chemkit::Real v1 = parameter(0);
    chemkit::Real v2 = parameter(1);
    chemkit::Real v3 = parameter(2);

    chemkit::Real phi = coordinates->torsionAngleRadians(a, b, c, d);

    // dE/dphi
    chemkit::Real de_dphi = (1.0/2.0) * (-v1 * sin(phi) + 2.0 * v2 * sin(2.0 * phi) - 3.0 * v3 * sin(3.0 * phi));

    boost::array<chemkit::Vector3, 4> gradient = coordinates->torsionAngleGradientRadians(a, b, c, d);

    gradient[0] *= de_dphi;
    gradient[1] *= de_dphi;
    gradient[2] *= de_dphi;
    gradient[3] *= de_dphi;

    return std::vector<chemkit::Vector3>(gradient.begin(), gradient.end());
}

// === OplsNonbondedCalculation ============================================ //
OplsNonbondedCalculation::OplsNonbondedCalculation(size_t a, size_t b)
    : OplsCalculation(VanDerWaals | Electrostatic, 2, 5)
{
    setAtom(0, a);
    setAtom(1, b);
}

bool OplsNonbondedCalculation::setup(const OplsParameters *parameters)
{
    int typeA = boost::lexical_cast<int>(atomType(0));
    int typeB = boost::lexical_cast<int>(atomType(1));

    const OplsVanDerWaalsParameters *pa = parameters->vanDerWaalsParameters(typeA);
    const OplsVanDerWaalsParameters *pb = parameters->vanDerWaalsParameters(typeB);
    if(!pa || !pb){
        return false;
    }

    chemkit::Real qa = parameters->partialCharge(typeA);
    chemkit::Real qb = parameters->partialCharge(typeB);
    chemkit::Real sigma = sqrt(pa->sigma * pb->sigma);
    chemkit::Real epsilon = sqrt(pa->epsilon * pb->epsilon);

    setParameter(0, qa);
    setParameter(1, qb);
    setParameter(2, sigma);
    setParameter(3, epsilon);

    // one-four scaling
    if(topology()->isOneFour(atom(0), atom(1))){
        setParameter(4, 0.5);
    }
    else{
        setParameter(4, 1.0);
    }

    return true;
}

chemkit::Real OplsNonbondedCalculation::energy(const chemkit::CartesianCoordinates *coordinates) const
{
    size_t a = atom(0);
    size_t b = atom(1);

    chemkit::Real qa = parameter(0);
    chemkit::Real qb = parameter(1);
    chemkit::Real e = 332.06; // vacuum permitivity
    chemkit::Real sigma = parameter(2);
    chemkit::Real epsilon = parameter(3);
    chemkit::Real scale = parameter(4);

    chemkit::Real r = coordinates->distance(a, b);

    return scale * ((qa * qb * e) / r + 4.0 * epsilon * (pow(sigma / r, 12) - pow(sigma / r, 6)));
}

std::vector<chemkit::Vector3> OplsNonbondedCalculation::gradient(const chemkit::CartesianCoordinates *coordinates) const
{
    std::vector<chemkit::Vector3> gradient(2);

    size_t a = atom(0);
    size_t b = atom(1);

    chemkit::Real qa = parameter(0);
    chemkit::Real qb = parameter(1);
    chemkit::Real e = 332.06; // vacuum permitivity
    chemkit::Real sigma = parameter(2);
    chemkit::Real epsilon = parameter(3);
    chemkit::Real scale = parameter(4);

    chemkit::Real r = coordinates->distance(a, b);
    chemkit::Real sr = sigma / r;

    // dE/dr
    chemkit::Real de_dr = scale * ((1.0 / pow(r, 3)) * (-qa * qb * e + -4.0 * epsilon * sigma * (12.0 * pow(sr, 11) - 6.0 * pow(sr, 5))));

    // dE/da
    chemkit::Vector3 de_da = (coordinates->position(a) - coordinates->position(b)) * de_dr;

    gradient[0] = de_da;
    gradient[1] = -de_da;

    return gradient;
}
