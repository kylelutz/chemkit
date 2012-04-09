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

#include "mmffcalculation.h"

#include <boost/lexical_cast.hpp>

#include <chemkit/topology.h>
#include <chemkit/constants.h>
#include <chemkit/forcefield.h>
#include <chemkit/cartesiancoordinates.h>

#include "mmffparameters.h"

// === MmffCalculation ===================================================== //
MmffCalculation::MmffCalculation(int type, int atomCount, int parameterCount)
    : ForceFieldCalculation(type, atomCount, parameterCount)
{
}

// === MmffBondStrechCalculation =========================================== //
MmffBondStrechCalculation::MmffBondStrechCalculation(size_t a, size_t b)
    : MmffCalculation(BondStrech, 2, 2)
{
    setAtom(0, a);
    setAtom(1, b);
}

bool MmffBondStrechCalculation::setup(const MmffParameters *parameters)
{
    const boost::shared_ptr<chemkit::Topology> &topology = this->topology();

    size_t a = atom(0);
    size_t b = atom(1);

    int typeA = boost::lexical_cast<int>(topology->type(a));
    int typeB = boost::lexical_cast<int>(topology->type(b));
    int bondType = topology->bondedInteractionType(a, b);

    const MmffBondStrechParameters *bondStrechParameters = parameters->bondStrechParameters(bondType, typeA, typeB);
    if(bondStrechParameters){
        setParameter(0, bondStrechParameters->kb);
        setParameter(1, bondStrechParameters->r0);
        return true;
    }

    return false;
}

chemkit::Real MmffBondStrechCalculation::energy(const chemkit::CartesianCoordinates *coordinates) const
{
    size_t a = atom(0);
    size_t b = atom(1);

    chemkit::Real kb = parameter(0);
    chemkit::Real r0 = parameter(1);

    chemkit::Real r = coordinates->distance(a, b);
    chemkit::Real dr = r - r0;
    chemkit::Real cs = -2.0; // cubic strech constant

    // equation 2
    return 143.9325 * (kb / 2) * (dr*dr) * (1 + cs * dr + ((7.0/12.0)*(cs*cs)) * (dr*dr));
}

std::vector<chemkit::Vector3> MmffBondStrechCalculation::gradient(const chemkit::CartesianCoordinates *coordinates) const
{
    size_t a = atom(0);
    size_t b = atom(1);

    chemkit::Real kb = parameter(0);
    chemkit::Real r0 = parameter(1);

    chemkit::Real r = coordinates->distance(a, b);
    chemkit::Real dr = r - r0;
    chemkit::Real cs = -2.0; // cubic strech constant

    // dE/dr
    chemkit::Real de_dr = 143.9325 * kb * dr * (1 + cs * dr + (7.0/12.0 * (cs*cs) * (dr*dr)) + 0.5 * dr * (cs + (14.0/12.0 * (cs*cs) * dr)));

    boost::array<chemkit::Vector3, 2> gradient = coordinates->distanceGradient(a, b);

    gradient[0] *= de_dr;
    gradient[1] *= de_dr;

    return std::vector<chemkit::Vector3>(gradient.begin(), gradient.end());
}

// === MmffAngleBendCalculation ============================================ //
MmffAngleBendCalculation::MmffAngleBendCalculation(size_t a, size_t b, size_t c)
    : MmffCalculation(AngleBend, 3, 2)
{
    setAtom(0, a);
    setAtom(1, b);
    setAtom(2, c);
}

bool MmffAngleBendCalculation::setup(const MmffParameters *parameters)
{
    const boost::shared_ptr<chemkit::Topology> &topology = this->topology();

    size_t a = atom(0);
    size_t b = atom(1);
    size_t c = atom(2);

    int typeA = boost::lexical_cast<int>(topology->type(a));
    int typeB = boost::lexical_cast<int>(topology->type(b));
    int typeC = boost::lexical_cast<int>(topology->type(c));
    int angleType = topology->angleInteractionType(a, b, c);

    const MmffAngleBendParameters *angleBendParameters =
        parameters->angleBendParameters(angleType, typeA, typeB, typeC);

    if(angleBendParameters){
        setParameter(0, angleBendParameters->ka);
        setParameter(1, angleBendParameters->theta0);
        return true;
    }

    return false;
}

chemkit::Real MmffAngleBendCalculation::energy(const chemkit::CartesianCoordinates *coordinates) const
{
    size_t a = atom(0);
    size_t b = atom(1);
    size_t c = atom(2);

    chemkit::Real ka = parameter(0);
    chemkit::Real t0 = parameter(1);

    chemkit::Real cb = -0.007; // cubic bend constant
    chemkit::Real t = coordinates->angle(a, b, c);
    chemkit::Real dt = t - t0;

    // equation 3
    return 0.043844 * (ka / 2.0) * pow(dt, 2) * (1 + cb * dt);
}

std::vector<chemkit::Vector3> MmffAngleBendCalculation::gradient(const chemkit::CartesianCoordinates *coordinates) const
{
    size_t a = atom(0);
    size_t b = atom(1);
    size_t c = atom(2);

    chemkit::Real ka = parameter(0);
    chemkit::Real t0 = parameter(1);

    chemkit::Real cb = -0.007; // cubic bend constant
    chemkit::Real t = coordinates->angle(a, b, c);
    chemkit::Real dt = t - t0;

    // dE/dt
    chemkit::Real de_dt = 0.043844 * ka * dt * (1 + cb * dt + 0.5 * cb * dt);

    boost::array<chemkit::Vector3, 3> gradient = coordinates->angleGradient(a, b, c);

    gradient[0] *= de_dt;
    gradient[1] *= de_dt;
    gradient[2] *= de_dt;

    return std::vector<chemkit::Vector3>(gradient.begin(), gradient.end());
}

// === MmffStrechBendCalculation =========================================== //
MmffStrechBendCalculation::MmffStrechBendCalculation(size_t a, size_t b, size_t c)
    : MmffCalculation(BondStrech | AngleBend, 3, 5)
{
    setAtom(0, a);
    setAtom(1, b);
    setAtom(2, c);
}

bool MmffStrechBendCalculation::setup(const MmffParameters *parameters)
{
    const boost::shared_ptr<chemkit::Topology> &topology = this->topology();

    size_t a = atom(0);
    size_t b = atom(1);
    size_t c = atom(2);

    int typeA = boost::lexical_cast<int>(topology->type(a));
    int typeB = boost::lexical_cast<int>(topology->type(b));
    int typeC = boost::lexical_cast<int>(topology->type(c));
    int bondTypeAB = topology->bondedInteractionType(a, b);
    int bondTypeBC = topology->bondedInteractionType(b, c);
    int angleType = topology->angleInteractionType(a, b, c);

    int strechBendType =
        parameters->calculateStrechBendType(bondTypeAB, bondTypeBC, angleType);

    bool parametersSwapped = false;
    const MmffStrechBendParameters *strechBendParameters =
        parameters->strechBendParameters(strechBendType, typeA, typeB, typeC);
    if(!strechBendParameters){
        strechBendType =
            parameters->calculateStrechBendType(bondTypeBC, bondTypeAB, angleType);
        strechBendParameters =
            parameters->strechBendParameters(strechBendType, typeC, typeB, typeA);

        if(strechBendParameters){
            parametersSwapped = true;
        }
        else{
            strechBendParameters = parameters->defaultStrechBendParameters(typeA, typeB, typeC);

            if(!strechBendParameters){
                strechBendParameters = parameters->defaultStrechBendParameters(typeC, typeB, typeA);
                parametersSwapped = true;
            }
        }
    }

    const MmffBondStrechParameters *bondStrechParameters_ab =
        parameters->bondStrechParameters(bondTypeAB, typeA, typeB);
    const MmffBondStrechParameters *bondStrechParameters_bc =
        parameters->bondStrechParameters(bondTypeBC, typeB, typeC);
    const MmffAngleBendParameters *angleBendParameters =
        parameters->angleBendParameters(angleType, typeA, typeB, typeC);
    if(strechBendParameters && bondStrechParameters_ab && bondStrechParameters_bc && angleBendParameters){
        if(parametersSwapped){
            setParameter(1, strechBendParameters->kba_ijk);
            setParameter(0, strechBendParameters->kba_kji);
        }
        else{
            setParameter(0, strechBendParameters->kba_ijk);
            setParameter(1, strechBendParameters->kba_kji);
        }

        setParameter(2, bondStrechParameters_ab->r0);
        setParameter(3, bondStrechParameters_bc->r0);
        setParameter(4, angleBendParameters->theta0);
        return true;
    }

    return false;
}

chemkit::Real MmffStrechBendCalculation::energy(const chemkit::CartesianCoordinates *coordinates) const
{
    size_t a = atom(0);
    size_t b = atom(1);
    size_t c = atom(2);

    chemkit::Real kba_ijk = parameter(0);
    chemkit::Real kba_kji = parameter(1);
    chemkit::Real r0_ab = parameter(2);
    chemkit::Real r0_bc = parameter(3);
    chemkit::Real t0 = parameter(4);

    chemkit::Real r_ab = coordinates->distance(a, b);
    chemkit::Real r_bc = coordinates->distance(b, c);
    chemkit::Real dr_ab = r_ab - r0_ab;
    chemkit::Real dr_bc = r_bc - r0_bc;
    chemkit::Real t = coordinates->angle(a, b, c);
    chemkit::Real dt = t - t0;

    // equation 5
    return 2.51210 * (kba_ijk * dr_ab + kba_kji * dr_bc) * dt;
}

std::vector<chemkit::Vector3> MmffStrechBendCalculation::gradient(const chemkit::CartesianCoordinates *coordinates) const
{
    size_t a = atom(0);
    size_t b = atom(1);
    size_t c = atom(2);

    chemkit::Real kba_ijk = parameter(0);
    chemkit::Real kba_kji = parameter(1);
    chemkit::Real r0_ab = parameter(2);
    chemkit::Real r0_bc = parameter(3);
    chemkit::Real t0 = parameter(4);

    chemkit::Real r_ab = coordinates->distance(a, b);
    chemkit::Real r_bc = coordinates->distance(b, c);
    chemkit::Real dr_ab = r_ab - r0_ab;
    chemkit::Real dr_bc = r_bc - r0_bc;
    chemkit::Real t = coordinates->angle(a, b, c);
    chemkit::Real dt = t - t0;

    std::vector<chemkit::Vector3> gradient(3);

    boost::array<chemkit::Vector3, 2> distanceGradientAB = coordinates->distanceGradient(a, b);
    boost::array<chemkit::Vector3, 2> distanceGradientBC = coordinates->distanceGradient(b, c);
    boost::array<chemkit::Vector3, 3> angleGradientABC = coordinates->angleGradient(a, b, c);

    gradient[0] = (distanceGradientAB[0] * kba_ijk * dt + angleGradientABC[0] * (kba_ijk * dr_ab + kba_kji * dr_bc)) * 2.51210;
    gradient[1] = ((distanceGradientAB[1] * kba_ijk + distanceGradientBC[0] * kba_kji) * dt + angleGradientABC[1] * (kba_ijk * dr_ab + kba_kji * dr_bc)) * 2.51210;
    gradient[2] = ((distanceGradientBC[1] * kba_kji) * dt + angleGradientABC[2] * (kba_ijk * dr_ab + kba_kji * dr_bc)) * 2.51210;

    return gradient;
}

// === MmffOutOfPlaneBendingCalculation ==================================== //
MmffOutOfPlaneBendingCalculation::MmffOutOfPlaneBendingCalculation(size_t a, size_t b, size_t c, size_t d)
    : MmffCalculation(Inversion, 4, 1)
{
    setAtom(0, a);
    setAtom(1, b);
    setAtom(2, c);
    setAtom(3, d);
}

bool MmffOutOfPlaneBendingCalculation::setup(const MmffParameters *parameters)
{
    const boost::shared_ptr<chemkit::Topology> &topology = this->topology();

    size_t a = atom(0);
    size_t b = atom(1);
    size_t c = atom(2);
    size_t d = atom(3);

    int typeA = boost::lexical_cast<int>(topology->type(a));
    int typeB = boost::lexical_cast<int>(topology->type(b));
    int typeC = boost::lexical_cast<int>(topology->type(c));
    int typeD = boost::lexical_cast<int>(topology->type(d));

    const MmffOutOfPlaneBendingParameters *outOfPlaneBendingParameters =
        parameters->outOfPlaneBendingParameters(typeA, typeB, typeC, typeD);
    if(!outOfPlaneBendingParameters){
        return false;
    }

    setParameter(0, outOfPlaneBendingParameters->koop);
    return true;
}

chemkit::Real MmffOutOfPlaneBendingCalculation::energy(const chemkit::CartesianCoordinates *coordinates) const
{
    size_t a = atom(0);
    size_t b = atom(1);
    size_t c = atom(2);
    size_t d = atom(3);

    chemkit::Real angle = coordinates->wilsonAngle(a, b, c, d);
    chemkit::Real koop = parameter(0);

    // equation 6
    return 0.043844 * (koop / 2.0) * (angle*angle);
}

std::vector<chemkit::Vector3> MmffOutOfPlaneBendingCalculation::gradient(const chemkit::CartesianCoordinates *coordinates) const
{
    size_t a = atom(0);
    size_t b = atom(1);
    size_t c = atom(2);
    size_t d = atom(3);

    chemkit::Real angle = coordinates->wilsonAngle(a, b, c, d);
    chemkit::Real koop = parameter(0);

    // dE/dw
    chemkit::Real de_dw = 0.043844 * koop * angle;

    boost::array<chemkit::Vector3, 4> gradient = coordinates->wilsonAngleGradient(a, b, c, d);

    gradient[0] *= de_dw;
    gradient[1] *= de_dw;
    gradient[2] *= de_dw;
    gradient[3] *= de_dw;

    return std::vector<chemkit::Vector3>(gradient.begin(), gradient.end());
}

// === MmffTorsionCalculation ============================================== //
MmffTorsionCalculation::MmffTorsionCalculation(size_t a, size_t b, size_t c, size_t d)
    : MmffCalculation(Torsion, 4, 3)
{
    setAtom(0, a);
    setAtom(1, b);
    setAtom(2, c);
    setAtom(3, d);
}

bool MmffTorsionCalculation::setup(const MmffParameters *parameters)
{
    const boost::shared_ptr<chemkit::Topology> &topology = this->topology();

    size_t a = atom(0);
    size_t b = atom(1);
    size_t c = atom(2);
    size_t d = atom(3);

    int typeA = boost::lexical_cast<int>(topology->type(a));
    int typeB = boost::lexical_cast<int>(topology->type(b));
    int typeC = boost::lexical_cast<int>(topology->type(c));
    int typeD = boost::lexical_cast<int>(topology->type(d));
    int torsionType = topology->torsionInteractionType(a, b, c, d);

    const MmffTorsionParameters *torsionParameters =
        parameters->torsionParameters(torsionType, typeA, typeB, typeC, typeD);
    if(!torsionParameters){
        return false;
    }

    setParameter(0, torsionParameters->V1);
    setParameter(1, torsionParameters->V2);
    setParameter(2, torsionParameters->V3);

    return true;
}

chemkit::Real MmffTorsionCalculation::energy(const chemkit::CartesianCoordinates *coordinates) const
{
    size_t a = atom(0);
    size_t b = atom(1);
    size_t c = atom(2);
    size_t d = atom(3);

    chemkit::Real angle = coordinates->torsionAngleRadians(a, b, c, d);
    chemkit::Real V1 = parameter(0);
    chemkit::Real V2 = parameter(1);
    chemkit::Real V3 = parameter(2);

    // equation 7
    return 0.5 * (V1 * (1.0 + cos(angle)) + V2 * (1.0 - cos(2.0 * angle)) + V3 * (1.0 + cos(3.0 * angle)));
}

std::vector<chemkit::Vector3> MmffTorsionCalculation::gradient(const chemkit::CartesianCoordinates *coordinates) const
{
    size_t a = atom(0);
    size_t b = atom(1);
    size_t c = atom(2);
    size_t d = atom(3);

    chemkit::Real phi = coordinates->torsionAngleRadians(a, b, c, d);
    chemkit::Real V1 = parameter(0);
    chemkit::Real V2 = parameter(1);
    chemkit::Real V3 = parameter(2);

    // dE/dphi
    chemkit::Real de_dphi = 0.5 * (-V1 * sin(phi) + 2 * V2 * sin(2 * phi) - 3 * V3 * sin(3 * phi));

    boost::array<chemkit::Vector3, 4> gradient = coordinates->torsionAngleGradientRadians(a, b, c, d);

    gradient[0] *= de_dphi;
    gradient[1] *= de_dphi;
    gradient[2] *= de_dphi;
    gradient[3] *= de_dphi;

    return std::vector<chemkit::Vector3>(gradient.begin(), gradient.end());
}

// === MmffVanDerWaalsCalculation ========================================== //
MmffVanDerWaalsCalculation::MmffVanDerWaalsCalculation(size_t a, size_t b)
    : MmffCalculation(VanDerWaals, 2, 2)
{
    setAtom(0, a);
    setAtom(1, b);
}

bool MmffVanDerWaalsCalculation::setup(const MmffParameters *parameters)
{
    const boost::shared_ptr<chemkit::Topology> &topology = this->topology();

    size_t a = atom(0);
    size_t b = atom(1);

    int typeA = boost::lexical_cast<int>(topology->type(a));
    int typeB = boost::lexical_cast<int>(topology->type(b));

    const MmffVanDerWaalsParameters *parametersA = parameters->vanDerWaalsParameters(typeA);
    const MmffVanDerWaalsParameters *parametersB = parameters->vanDerWaalsParameters(typeB);
    if(!parametersA || !parametersB){
        return false;
    }

    chemkit::Real N_a = parametersA->N;
    chemkit::Real N_b = parametersB->N;
    chemkit::Real A_a = parametersA->A;
    chemkit::Real A_b = parametersB->A;
    chemkit::Real G_a = parametersA->G;
    chemkit::Real G_b = parametersB->G;
    chemkit::Real alpha_a = parametersA->alpha;
    chemkit::Real alpha_b = parametersB->alpha;
    char DA_a = parametersA->DA;
    char DA_b = parametersB->DA;

    // equation 9
    chemkit::Real rs_aa = A_a * pow(alpha_a, (1.0/4.0));
    chemkit::Real rs_bb = A_b * pow(alpha_b, (1.0/4.0));

    // equation 11
    chemkit::Real gamma = (rs_aa - rs_bb) / (rs_aa + rs_bb);

    // equation 10
    chemkit::Real rs;

    if(DA_a == 'D' || (DA_b == 'D')){
        rs = 0.5 * (rs_aa + rs_bb);
    }
    else{
        rs = 0.5 * (rs_aa + rs_bb) * (1.0 + 0.2 * (1.0 - exp(-12.0 * gamma * gamma)));
    }

    // equation 12
    chemkit::Real eps = ((181.16 * G_a * G_b * alpha_a * alpha_b) / (sqrt(alpha_a / N_a) + sqrt(alpha_b / N_b))) * pow(rs, -6.0);

    if((DA_a == 'D' && DA_b == 'A') || (DA_a == 'A' && DA_b == 'D')){
        rs *= 0.8;
        eps *= 0.5;
    }

    setParameter(0, rs);
    setParameter(1, eps);

    return true;
}

chemkit::Real MmffVanDerWaalsCalculation::energy(const chemkit::CartesianCoordinates *coordinates) const
{
    size_t a = atom(0);
    size_t b = atom(1);

    chemkit::Real rs = parameter(0);
    chemkit::Real eps = parameter(1);
    chemkit::Real r = coordinates->distance(a, b);

    // equation 8
    return eps * pow(((1.07 * rs) / (r + 0.07 * rs)), 7) * (((1.12 * pow(rs, 7)) / (pow(r, 7) + 0.12 * pow(rs, 7))) - 2);
}

std::vector<chemkit::Vector3> MmffVanDerWaalsCalculation::gradient(const chemkit::CartesianCoordinates *coordinates) const
{
    size_t a = atom(0);
    size_t b = atom(1);

    chemkit::Real rs = parameter(0);
    chemkit::Real eps = parameter(1);
    chemkit::Real r = coordinates->distance(a, b);

    // dE/dr
    chemkit::Real de_dr = 7 * eps * pow(1.07 * rs / (r + 0.07 * rs), 6) *
                           ((-1.07 * rs / pow(r + 0.07 * rs, 2)) * (1.12 * pow(rs, 7) / (pow(r, 7) + 0.12 * pow(rs, 7)) - 2) +
                           (-1.12 * pow(rs, 7) * pow(r, 6) / pow(pow(r, 7) + 0.12 * pow(rs, 7), 2)) * (1.07 * rs / (r + 0.07 * rs)));

    boost::array<chemkit::Vector3, 2> gradient = coordinates->distanceGradient(a, b);

    gradient[0] *= de_dr;
    gradient[1] *= de_dr;

    return std::vector<chemkit::Vector3>(gradient.begin(), gradient.end());
}

// === MmffElectrostaticCalculation ======================================== //
MmffElectrostaticCalculation::MmffElectrostaticCalculation(size_t a, size_t b)
    : MmffCalculation(Electrostatic, 2, 3)
{
    setAtom(0, a);
    setAtom(1, b);
}

bool MmffElectrostaticCalculation::setup(const MmffParameters *parameters)
{
    CHEMKIT_UNUSED(parameters);

    const boost::shared_ptr<chemkit::Topology> &topology = this->topology();

    size_t a = atom(0);
    size_t b = atom(1);

    chemkit::Real oneFourScaling;

    if(topology->isOneFour(a, b)){
        oneFourScaling = 0.75;
    }
    else{
        oneFourScaling = 1.0;
    }

    setParameter(0, topology->charge(a));
    setParameter(1, topology->charge(b));
    setParameter(2, oneFourScaling);

    return true;
}

chemkit::Real MmffElectrostaticCalculation::energy(const chemkit::CartesianCoordinates *coordinates) const
{
    size_t a = atom(0);
    size_t b = atom(1);

    chemkit::Real qa = parameter(0);
    chemkit::Real qb = parameter(1);
    chemkit::Real oneFourScaling = parameter(2);

    chemkit::Real r = coordinates->distance(a, b);
    chemkit::Real e = 1.0; // dielectric constant
    chemkit::Real d = 0.05; // electrostatic buffering constant

    // equation 13
    return ((332.0716 * qa * qb) / (e * (r + d))) * oneFourScaling;
}

std::vector<chemkit::Vector3> MmffElectrostaticCalculation::gradient(const chemkit::CartesianCoordinates *coordinates) const
{
    size_t a = atom(0);
    size_t b = atom(1);

    chemkit::Real qa = parameter(0);
    chemkit::Real qb = parameter(1);
    chemkit::Real oneFourScaling = parameter(2);

    chemkit::Real r = coordinates->distance(a, b);
    chemkit::Real e = 1.0; // dielectric constant
    chemkit::Real d = 0.05; // electrostatic buffering constant

    chemkit::Real de_dr = 332.0716 * qa * qb * oneFourScaling * (-1.0 / (e * pow(r + d, 2)));

    boost::array<chemkit::Vector3, 2> gradient = coordinates->distanceGradient(a, b);

    gradient[0] *= de_dr;
    gradient[1] *= de_dr;

    return std::vector<chemkit::Vector3>(gradient.begin(), gradient.end());
}
