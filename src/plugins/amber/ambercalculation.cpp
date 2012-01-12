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

#include "ambercalculation.h"

#include "amberparameters.h"

#include <chemkit/constants.h>
#include <chemkit/forcefieldatom.h>

// === AmberCalculation ==================================================== //
AmberCalculation::AmberCalculation(int type, int atomCount, int parameterCount)
    : ForceFieldCalculation(type, atomCount, parameterCount)
{
}

// === AmberBondCalculation ================================================ //
AmberBondCalculation::AmberBondCalculation(const chemkit::ForceFieldAtom *a,
                                           const chemkit::ForceFieldAtom *b)
    : AmberCalculation(BondStrech, 2, 2)
{
    setAtom(0, a);
    setAtom(1, b);
}

bool AmberBondCalculation::setup(const AmberParameters *parameters)
{
    const chemkit::ForceFieldAtom *a = atom(0);
    const chemkit::ForceFieldAtom *b = atom(1);

    const AmberBondParameters *bondParameters = parameters->bondParameters(a, b);
    if(!bondParameters){
        return false;
    }

    setParameter(0, bondParameters->kb);
    setParameter(1, bondParameters->r0);

    return true;
}

chemkit::Real AmberBondCalculation::energy() const
{
    const chemkit::ForceFieldAtom *a = atom(0);
    const chemkit::ForceFieldAtom *b = atom(1);

    chemkit::Real kb = parameter(0);
    chemkit::Real r0 = parameter(1);
    chemkit::Real r = distance(a, b);
    chemkit::Real dr = r - r0;

    return kb * (dr*dr);
}

std::vector<chemkit::Vector3> AmberBondCalculation::gradient() const
{
    const chemkit::ForceFieldAtom *a = atom(0);
    const chemkit::ForceFieldAtom *b = atom(1);

    chemkit::Real kb = parameter(0);
    chemkit::Real r0 = parameter(1);
    chemkit::Real r = distance(a, b);

    // dE/dr
    chemkit::Real de_dr = 2.0 * kb * (r - r0);

    boost::array<chemkit::Vector3, 2> gradient = distanceGradient(a, b);

    gradient[0] *= de_dr;
    gradient[1] *= de_dr;

    return std::vector<chemkit::Vector3>(gradient.begin(), gradient.end());
}

// === AmberAngleCalculation =============================================== //
AmberAngleCalculation::AmberAngleCalculation(const chemkit::ForceFieldAtom *a,
                                             const chemkit::ForceFieldAtom *b,
                                             const chemkit::ForceFieldAtom *c)
    : AmberCalculation(AngleBend, 3, 2)
{
    setAtom(0, a);
    setAtom(1, b);
    setAtom(2, c);
}

bool AmberAngleCalculation::setup(const AmberParameters *parameters)
{
    const chemkit::ForceFieldAtom *a = atom(0);
    const chemkit::ForceFieldAtom *b = atom(1);
    const chemkit::ForceFieldAtom *c = atom(2);

    const AmberAngleParameters *angleParameters = parameters->angleParameters(a, b, c);
    if(!angleParameters){
        return false;
    }

    setParameter(0, angleParameters->ka);
    setParameter(1, angleParameters->theta0);

    return true;
}

chemkit::Real AmberAngleCalculation::energy() const
{
    const chemkit::ForceFieldAtom *a = atom(0);
    const chemkit::ForceFieldAtom *b = atom(1);
    const chemkit::ForceFieldAtom *c = atom(2);

    chemkit::Real ka = parameter(0);
    chemkit::Real theta0 = parameter(1);
    chemkit::Real theta = bondAngle(a, b, c);
    chemkit::Real dt = theta - theta0;

    return ka * (dt*dt);
}

std::vector<chemkit::Vector3> AmberAngleCalculation::gradient() const
{
    const chemkit::ForceFieldAtom *a = atom(0);
    const chemkit::ForceFieldAtom *b = atom(1);
    const chemkit::ForceFieldAtom *c = atom(2);

    chemkit::Real ka = parameter(0);
    chemkit::Real theta0 = parameter(1);
    chemkit::Real theta = bondAngle(a, b, c);

    // dE/dtheta
    chemkit::Real de_dtheta = 2.0 * ka * (theta - theta0);

    boost::array<chemkit::Vector3, 3> gradient = bondAngleGradient(a, b, c);

    gradient[0] *= de_dtheta;
    gradient[1] *= de_dtheta;
    gradient[2] *= de_dtheta;

    return std::vector<chemkit::Vector3>(gradient.begin(), gradient.end());
}

// === AmberTorsionCalculation ============================================= //
AmberTorsionCalculation::AmberTorsionCalculation(const chemkit::ForceFieldAtom *a,
                                                 const chemkit::ForceFieldAtom *b,
                                                 const chemkit::ForceFieldAtom *c,
                                                 const chemkit::ForceFieldAtom *d)
    : AmberCalculation(Torsion, 4, 8)
{
    setAtom(0, a);
    setAtom(1, b);
    setAtom(2, c);
    setAtom(3, d);
}

bool AmberTorsionCalculation::setup(const AmberParameters *parameters)
{
    const chemkit::ForceFieldAtom *a = atom(0);
    const chemkit::ForceFieldAtom *b = atom(1);
    const chemkit::ForceFieldAtom *c = atom(2);
    const chemkit::ForceFieldAtom *d = atom(3);

    const AmberTorsionParameters *torsionParameters = parameters->torsionParameters(a, b, c, d);
    if(!torsionParameters){
        return false;
    }

    setParameter(0, torsionParameters->V1);
    setParameter(1, torsionParameters->V2);
    setParameter(2, torsionParameters->V3);
    setParameter(3, torsionParameters->V4);
    setParameter(4, torsionParameters->gamma1);
    setParameter(5, torsionParameters->gamma2);
    setParameter(6, torsionParameters->gamma3);
    setParameter(7, torsionParameters->gamma4);

    return true;
}

chemkit::Real AmberTorsionCalculation::energy() const
{
    const chemkit::ForceFieldAtom *a = atom(0);
    const chemkit::ForceFieldAtom *b = atom(1);
    const chemkit::ForceFieldAtom *c = atom(2);
    const chemkit::ForceFieldAtom *d = atom(3);

    chemkit::Real V1 = parameter(0);
    chemkit::Real V2 = parameter(1);
    chemkit::Real V3 = parameter(2);
    chemkit::Real V4 = parameter(3);
    chemkit::Real gamma1 = parameter(4);
    chemkit::Real gamma2 = parameter(5);
    chemkit::Real gamma3 = parameter(6);
    chemkit::Real gamma4 = parameter(7);

    chemkit::Real angle = torsionAngle(a, b, c, d);

    chemkit::Real energy = 0;
    energy += V1 * (1.0 + cos((1.0 * angle - gamma1) * chemkit::constants::DegreesToRadians));
    energy += V2 * (1.0 + cos((2.0 * angle - gamma2) * chemkit::constants::DegreesToRadians));
    energy += V3 * (1.0 + cos((3.0 * angle - gamma3) * chemkit::constants::DegreesToRadians));
    energy += V4 * (1.0 + cos((4.0 * angle - gamma4) * chemkit::constants::DegreesToRadians));

    return energy;
}

std::vector<chemkit::Vector3> AmberTorsionCalculation::gradient() const
{
    const chemkit::ForceFieldAtom *a = atom(0);
    const chemkit::ForceFieldAtom *b = atom(1);
    const chemkit::ForceFieldAtom *c = atom(2);
    const chemkit::ForceFieldAtom *d = atom(3);

    chemkit::Real V1 = parameter(0);
    chemkit::Real V2 = parameter(1);
    chemkit::Real V3 = parameter(2);
    chemkit::Real V4 = parameter(3);
    chemkit::Real gamma1 = parameter(4);
    chemkit::Real gamma2 = parameter(5);
    chemkit::Real gamma3 = parameter(6);
    chemkit::Real gamma4 = parameter(7);

    chemkit::Real phi = torsionAngle(a, b, c, d);

    // dE/dphi
    chemkit::Real de_dphi = 0;
    de_dphi += V1 * (-sin((1.0 * phi - gamma1) * chemkit::constants::DegreesToRadians) * 1.0);
    de_dphi += V2 * (-sin((2.0 * phi - gamma2) * chemkit::constants::DegreesToRadians) * 2.0);
    de_dphi += V3 * (-sin((3.0 * phi - gamma3) * chemkit::constants::DegreesToRadians) * 3.0);
    de_dphi += V4 * (-sin((4.0 * phi - gamma4) * chemkit::constants::DegreesToRadians) * 4.0);
    de_dphi *= chemkit::constants::DegreesToRadians;

    boost::array<chemkit::Vector3, 4> gradient = torsionAngleGradient(a, b, c, d);

    gradient[0] *= de_dphi;
    gradient[1] *= de_dphi;
    gradient[2] *= de_dphi;
    gradient[3] *= de_dphi;

    return std::vector<chemkit::Vector3>(gradient.begin(), gradient.end());
}

// === AmberNonbondedCalculation =========================================== //
AmberNonbondedCalculation::AmberNonbondedCalculation(const chemkit::ForceFieldAtom *a,
                                                     const chemkit::ForceFieldAtom *b)
    : AmberCalculation(VanDerWaals | Electrostatic, 2, 2)
{
    setAtom(0, a);
    setAtom(1, b);
}

bool AmberNonbondedCalculation::setup(const AmberParameters *parameters)
{
    const chemkit::ForceFieldAtom *a = atom(0);
    const chemkit::ForceFieldAtom *b = atom(1);

    const struct AmberNonbondedParameters *parametersA = parameters->nonbondedParameters(a);
    const struct AmberNonbondedParameters *parametersB = parameters->nonbondedParameters(b);
    if(!parametersA || !parametersB){
        return false;
    }

    chemkit::Real epsilon = parametersA->wellDepth + parametersB->wellDepth;
    chemkit::Real sigma = parametersA->vanDerWaalsRadius + parametersB->vanDerWaalsRadius;

    setParameter(0, epsilon);
    setParameter(1, sigma);

    return true;
}

chemkit::Real AmberNonbondedCalculation::energy() const
{
    const chemkit::ForceFieldAtom *a = atom(0);
    const chemkit::ForceFieldAtom *b = atom(1);

    chemkit::Real epsilon = parameter(0);
    chemkit::Real sigma = parameter(1);
    chemkit::Real r = distance(a, b);
    chemkit::Real e0 = 1;

    chemkit::Real vanDerWaalsTerm = epsilon * (pow(sigma/r, 12) - 2 * pow(sigma/r, 6));
    chemkit::Real electrostaticTerm = (a->charge() * b->charge()) / (4.0 * chemkit::constants::Pi * e0 * r);

    return vanDerWaalsTerm + electrostaticTerm;
}

std::vector<chemkit::Vector3> AmberNonbondedCalculation::gradient() const
{
    const chemkit::ForceFieldAtom *a = atom(0);
    const chemkit::ForceFieldAtom *b = atom(1);

    chemkit::Real epsilon = parameter(0);
    chemkit::Real sigma = parameter(1);
    chemkit::Real qa = a->charge();
    chemkit::Real qb = b->charge();
    chemkit::Real e0 = 1;
    chemkit::Real pi = chemkit::constants::Pi;

    chemkit::Real r = distance(a, b);
    chemkit::Real sr = sigma / r;

    // dE/dr
    chemkit::Real de_dr = (-12 * epsilon * sigma / pow(r, 2) * (pow(sr, 11) - pow(sr, 5))) - ((qa * qb) / (4.0 * pi * e0 * pow(r, 2)));

    boost::array<chemkit::Vector3, 2> gradient = distanceGradient(a, b);

    gradient[0] *= de_dr;
    gradient[1] *= de_dr;

    return std::vector<chemkit::Vector3>(gradient.begin(), gradient.end());
}
