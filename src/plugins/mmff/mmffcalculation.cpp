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

#include "mmffcalculation.h"

#include <chemkit/atom.h>
#include <chemkit/constants.h>
#include <chemkit/forcefield.h>

#include "mmffatom.h"
#include "mmffparameters.h"

// === MmffCalculation ===================================================== //
MmffCalculation::MmffCalculation(int type, int atomCount, int parameterCount)
    : ForceFieldCalculation(type, atomCount, parameterCount)
{
}

const MmffAtom* MmffCalculation::atom(int index) const
{
    return static_cast<const MmffAtom *>(ForceFieldCalculation::atom(index));
}

// === MmffBondStrechCalculation =========================================== //
MmffBondStrechCalculation::MmffBondStrechCalculation(const MmffAtom *a,
                                                     const MmffAtom *b)
    : MmffCalculation(BondStrech, 2, 2)
{
    setAtom(0, a);
    setAtom(1, b);
}

bool MmffBondStrechCalculation::setup(const MmffParameters *parameters)
{
    const MmffAtom *a = atom(0);
    const MmffAtom *b = atom(1);

    const MmffBondStrechParameters *bondStrechParameters = parameters->bondStrechParameters(a, b);
    if(bondStrechParameters){
        setParameter(0, bondStrechParameters->kb);
        setParameter(1, bondStrechParameters->r0);
        return true;
    }

    return false;
}

chemkit::Float MmffBondStrechCalculation::energy() const
{
    const MmffAtom *a = atom(0);
    const MmffAtom *b = atom(1);

    chemkit::Float kb = parameter(0);
    chemkit::Float r0 = parameter(1);

    chemkit::Float r = distance(a, b);
    chemkit::Float dr = r - r0;
    chemkit::Float cs = -2.0; // cubic strech constant

    // equation 2
    return 143.9325 * (kb / 2) * (dr*dr) * (1 + cs * dr + ((7.0/12.0)*(cs*cs)) * (dr*dr));
}

QVector<chemkit::Vector3> MmffBondStrechCalculation::gradient() const
{
    const MmffAtom *a = atom(0);
    const MmffAtom *b = atom(1);

    chemkit::Float kb = parameter(0);
    chemkit::Float r0 = parameter(1);

    chemkit::Float r = distance(a, b);
    chemkit::Float dr = r - r0;
    chemkit::Float cs = -2.0; // cubic strech constant

    // dE/dr
    chemkit::Float de_dr = 143.9325 * kb * dr * (1 + cs * dr + (7.0/12.0 * (cs*cs) * (dr*dr)) + 0.5 * dr * (cs + (14.0/12.0 * (cs*cs) * dr)));

    QVector<chemkit::Vector3> gradient = distanceGradient(a, b);

    gradient[0] *= de_dr;
    gradient[1] *= de_dr;

    return gradient;
}

// === MmffAngleBendCalculation ============================================ //
MmffAngleBendCalculation::MmffAngleBendCalculation(const MmffAtom *a,
                                                   const MmffAtom *b,
                                                   const MmffAtom *c)
    : MmffCalculation(AngleBend, 3, 2)
{
    setAtom(0, a);
    setAtom(1, b);
    setAtom(2, c);
}

bool MmffAngleBendCalculation::setup(const MmffParameters *parameters)
{
    const MmffAtom *a = atom(0);
    const MmffAtom *b = atom(1);
    const MmffAtom *c = atom(2);

    const MmffAngleBendParameters *angleBendParameters = parameters->angleBendParameters(a, b, c);
    if(angleBendParameters){
        setParameter(0, angleBendParameters->ka);
        setParameter(1, angleBendParameters->theta0);
        return true;
    }

    return false;
}

chemkit::Float MmffAngleBendCalculation::energy() const
{
    const MmffAtom *a = atom(0);
    const MmffAtom *b = atom(1);
    const MmffAtom *c = atom(2);

    chemkit::Float ka = parameter(0);
    chemkit::Float t0 = parameter(1);

    chemkit::Float cb = -0.007; // cubic bend constant
    chemkit::Float t = bondAngle(a, b, c);
    chemkit::Float dt = t - t0;

    // equation 3
    return 0.043844 * (ka / 2.0) * pow(dt, 2) * (1 + cb * dt);
}

QVector<chemkit::Vector3> MmffAngleBendCalculation::gradient() const
{
    const MmffAtom *a = atom(0);
    const MmffAtom *b = atom(1);
    const MmffAtom *c = atom(2);

    chemkit::Float ka = parameter(0);
    chemkit::Float t0 = parameter(1);

    chemkit::Float cb = -0.007; // cubic bend constant
    chemkit::Float t = bondAngle(a, b, c);
    chemkit::Float dt = t - t0;

    // dE/dt
    chemkit::Float de_dt = 0.043844 * ka * dt * (1 + cb * dt + 0.5 * cb * dt);

    QVector<chemkit::Vector3> gradient = bondAngleGradient(a, b, c);

    gradient[0] *= de_dt;
    gradient[1] *= de_dt;
    gradient[2] *= de_dt;

    return gradient;
}

// === MmffStrechBendCalculation =========================================== //
MmffStrechBendCalculation::MmffStrechBendCalculation(const MmffAtom *a,
                                                     const MmffAtom *b,
                                                     const MmffAtom *c)
    : MmffCalculation(BondStrech | AngleBend, 3, 5)
{
    setAtom(0, a);
    setAtom(1, b);
    setAtom(2, c);
}

bool MmffStrechBendCalculation::setup(const MmffParameters *parameters)
{
    const MmffAtom *a = atom(0);
    const MmffAtom *b = atom(1);
    const MmffAtom *c = atom(2);

    bool parametersSwapped = false;
    const MmffStrechBendParameters *strechBendParameters = parameters->strechBendParameters(a, b, c);
    if(!strechBendParameters){
        strechBendParameters = parameters->strechBendParameters(c, b, a);

        if(strechBendParameters){
            parametersSwapped = true;
        }
        else{
            strechBendParameters = parameters->defaultStrechBendParameters(a, b, c);

            if(!strechBendParameters){
                strechBendParameters = parameters->defaultStrechBendParameters(c, b, a);
                parametersSwapped = true;
            }
        }
    }

    const MmffBondStrechParameters *bondStrechParameters_ab = parameters->bondStrechParameters(a, b);
    const MmffBondStrechParameters *bondStrechParameters_bc = parameters->bondStrechParameters(b, c);
    const MmffAngleBendParameters *angleBendParameters = parameters->angleBendParameters(a, b, c);
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

chemkit::Float MmffStrechBendCalculation::energy() const
{
    const MmffAtom *a = atom(0);
    const MmffAtom *b = atom(1);
    const MmffAtom *c = atom(2);

    chemkit::Float kba_ijk = parameter(0);
    chemkit::Float kba_kji = parameter(1);
    chemkit::Float r0_ab = parameter(2);
    chemkit::Float r0_bc = parameter(3);
    chemkit::Float t0 = parameter(4);

    chemkit::Float r_ab = distance(a, b);
    chemkit::Float r_bc = distance(b, c);
    chemkit::Float dr_ab = r_ab - r0_ab;
    chemkit::Float dr_bc = r_bc - r0_bc;
    chemkit::Float t = bondAngle(a, b, c);
    chemkit::Float dt = t - t0;

    // equation 5
    return 2.51210 * (kba_ijk * dr_ab + kba_kji * dr_bc) * dt;
}

QVector<chemkit::Vector3> MmffStrechBendCalculation::gradient() const
{
    const MmffAtom *a = atom(0);
    const MmffAtom *b = atom(1);
    const MmffAtom *c = atom(2);

    chemkit::Float kba_ijk = parameter(0);
    chemkit::Float kba_kji = parameter(1);
    chemkit::Float r0_ab = parameter(2);
    chemkit::Float r0_bc = parameter(3);
    chemkit::Float t0 = parameter(4);

    chemkit::Float r_ab = distance(a, b);
    chemkit::Float r_bc = distance(b, c);
    chemkit::Float dr_ab = r_ab - r0_ab;
    chemkit::Float dr_bc = r_bc - r0_bc;
    chemkit::Float t = bondAngle(a, b, c);
    chemkit::Float dt = t - t0;

    QVector<chemkit::Vector3> gradient(3);

    QVector<chemkit::Vector3> distanceGradientAB = distanceGradient(a, b);
    QVector<chemkit::Vector3> distanceGradientBC = distanceGradient(b, c);
    QVector<chemkit::Vector3> bondAngleGradientABC = bondAngleGradient(a, b, c);

    gradient[0] = (distanceGradientAB[0] * kba_ijk * dt + bondAngleGradientABC[0] * (kba_ijk * dr_ab + kba_kji * dr_bc)) * 2.51210;
    gradient[1] = ((distanceGradientAB[1] * kba_ijk + distanceGradientBC[0] * kba_kji) * dt + bondAngleGradientABC[1] * (kba_ijk * dr_ab + kba_kji * dr_bc)) * 2.51210;
    gradient[2] = ((distanceGradientBC[1] * kba_kji) * dt + bondAngleGradientABC[2] * (kba_ijk * dr_ab + kba_kji * dr_bc)) * 2.51210;

    return gradient;
}

// === MmffOutOfPlaneBendingCalculation ==================================== //
MmffOutOfPlaneBendingCalculation::MmffOutOfPlaneBendingCalculation(const MmffAtom *a,
                                                                   const MmffAtom *b,
                                                                   const MmffAtom *c,
                                                                   const MmffAtom *d)
    : MmffCalculation(Inversion, 4, 1)
{
    setAtom(0, a);
    setAtom(1, b);
    setAtom(2, c);
    setAtom(3, d);
}

bool MmffOutOfPlaneBendingCalculation::setup(const MmffParameters *parameters)
{
    const MmffAtom *a = atom(0);
    const MmffAtom *b = atom(1);
    const MmffAtom *c = atom(2);
    const MmffAtom *d = atom(3);

    const MmffOutOfPlaneBendingParameters *outOfPlaneBendingParameters = parameters->outOfPlaneBendingParameters(a, b, c, d);
    if(!outOfPlaneBendingParameters){
        return false;
    }

    setParameter(0, outOfPlaneBendingParameters->koop);
    return true;
}

chemkit::Float MmffOutOfPlaneBendingCalculation::energy() const
{
    const MmffAtom *a = atom(0);
    const MmffAtom *b = atom(1);
    const MmffAtom *c = atom(2);
    const MmffAtom *d = atom(3);

    chemkit::Float angle = wilsonAngle(a, b, c, d);
    chemkit::Float koop = parameter(0);

    // equation 6
    return 0.043844 * (koop / 2.0) * (angle*angle);
}

QVector<chemkit::Vector3> MmffOutOfPlaneBendingCalculation::gradient() const
{
    const MmffAtom *a = atom(0);
    const MmffAtom *b = atom(1);
    const MmffAtom *c = atom(2);
    const MmffAtom *d = atom(3);

    chemkit::Float angle = wilsonAngle(a, b, c, d);
    chemkit::Float koop = parameter(0);

    // dE/dw
    chemkit::Float de_dw = 0.043844 * koop * angle;

    QVector<chemkit::Vector3> gradient = wilsonAngleGradient(a, b, c, d);

    gradient[0] *= de_dw;
    gradient[1] *= de_dw;
    gradient[2] *= de_dw;
    gradient[3] *= de_dw;

    return gradient;
}

// === MmffTorsionCalculation ============================================== //
MmffTorsionCalculation::MmffTorsionCalculation(const MmffAtom *a,
                                               const MmffAtom *b,
                                               const MmffAtom *c,
                                               const MmffAtom *d)
    : MmffCalculation(Torsion, 4, 3)
{
    setAtom(0, a);
    setAtom(1, b);
    setAtom(2, c);
    setAtom(3, d);
}

bool MmffTorsionCalculation::setup(const MmffParameters *parameters)
{
    const MmffAtom *a = atom(0);
    const MmffAtom *b = atom(1);
    const MmffAtom *c = atom(2);
    const MmffAtom *d = atom(3);

    const MmffTorsionParameters *torsionParameters = parameters->torsionParameters(a, b, c, d);
    if(!torsionParameters){
        return false;
    }

    setParameter(0, torsionParameters->V1);
    setParameter(1, torsionParameters->V2);
    setParameter(2, torsionParameters->V3);

    return true;
}

chemkit::Float MmffTorsionCalculation::energy() const
{
    const MmffAtom *a = atom(0);
    const MmffAtom *b = atom(1);
    const MmffAtom *c = atom(2);
    const MmffAtom *d = atom(3);

    chemkit::Float angle = torsionAngleRadians(a, b, c, d);
    chemkit::Float V1 = parameter(0);
    chemkit::Float V2 = parameter(1);
    chemkit::Float V3 = parameter(2);

    // equation 7
    return 0.5 * (V1 * (1.0 + cos(angle)) + V2 * (1.0 - cos(2.0 * angle)) + V3 * (1.0 + cos(3.0 * angle)));
}

QVector<chemkit::Vector3> MmffTorsionCalculation::gradient() const
{
    const MmffAtom *a = atom(0);
    const MmffAtom *b = atom(1);
    const MmffAtom *c = atom(2);
    const MmffAtom *d = atom(3);

    chemkit::Float phi = torsionAngleRadians(a, b, c, d);
    chemkit::Float V1 = parameter(0);
    chemkit::Float V2 = parameter(1);
    chemkit::Float V3 = parameter(2);

    // dE/dphi
    chemkit::Float de_dphi = 0.5 * (-V1 * sin(phi) + 2 * V2 * sin(2 * phi) - 3 * V3 * sin(3 * phi));

    QVector<chemkit::Vector3> gradient = torsionAngleGradientRadians(a, b, c, d);

    gradient[0] *= de_dphi;
    gradient[1] *= de_dphi;
    gradient[2] *= de_dphi;
    gradient[3] *= de_dphi;

    return gradient;
}

// === MmffVanDerWaalsCalculation ========================================== //
MmffVanDerWaalsCalculation::MmffVanDerWaalsCalculation(const MmffAtom *a,
                                                       const MmffAtom *b)
    : MmffCalculation(VanDerWaals, 2, 2)
{
    setAtom(0, a);
    setAtom(1, b);
}

bool MmffVanDerWaalsCalculation::setup(const MmffParameters *parameters)
{
    const MmffAtom *a = atom(0);
    const MmffAtom *b = atom(1);

    const MmffVanDerWaalsParameters *parametersA = parameters->vanDerWaalsParameters(a);
    const MmffVanDerWaalsParameters *parametersB = parameters->vanDerWaalsParameters(b);
    if(!parametersA || !parametersB){
        return false;
    }

    chemkit::Float N_a = parametersA->N;
    chemkit::Float N_b = parametersB->N;
    chemkit::Float A_a = parametersA->A;
    chemkit::Float A_b = parametersB->A;
    chemkit::Float G_a = parametersA->G;
    chemkit::Float G_b = parametersB->G;
    chemkit::Float alpha_a = parametersA->alpha;
    chemkit::Float alpha_b = parametersB->alpha;
    char DA_a = parametersA->DA;
    char DA_b = parametersB->DA;

    // equation 9
    chemkit::Float rs_aa = A_a * pow(alpha_a, (1.0/4.0));
    chemkit::Float rs_bb = A_b * pow(alpha_b, (1.0/4.0));

    // equation 11
    chemkit::Float gamma = (rs_aa - rs_bb) / (rs_aa + rs_bb);

    // equation 10
    chemkit::Float rs;

    if(DA_a == 'D' || (DA_b == 'D')){
        rs = 0.5 * (rs_aa + rs_bb);
    }
    else{
        rs = 0.5 * (rs_aa + rs_bb) * (1.0 + 0.2 * (1.0 - exp(-12.0 * gamma * gamma)));
    }

    // equation 12
    chemkit::Float eps = ((181.16 * G_a * G_b * alpha_a * alpha_b) / (sqrt(alpha_a / N_a) + sqrt(alpha_b / N_b))) * pow(rs, -6.0);

    if((DA_a == 'D' && DA_b == 'A') || (DA_a == 'A' && DA_b == 'D')){
        rs *= 0.8;
        eps *= 0.5;
    }

    setParameter(0, rs);
    setParameter(1, eps);

    return true;
}

chemkit::Float MmffVanDerWaalsCalculation::energy() const
{
    const MmffAtom *a = atom(0);
    const MmffAtom *b = atom(1);

    chemkit::Float rs = parameter(0);
    chemkit::Float eps = parameter(1);
    chemkit::Float r = distance(a, b);

    // equation 8
    return eps * pow(((1.07 * rs) / (r + 0.07 * rs)), 7) * (((1.12 * pow(rs, 7)) / (pow(r, 7) + 0.12 * pow(rs, 7))) - 2);
}

QVector<chemkit::Vector3> MmffVanDerWaalsCalculation::gradient() const
{
    const MmffAtom *a = atom(0);
    const MmffAtom *b = atom(1);

    chemkit::Float rs = parameter(0);
    chemkit::Float eps = parameter(1);
    chemkit::Float r = distance(a, b);

    // dE/dr
    chemkit::Float de_dr = 7 * eps * pow(1.07 * rs / (r + 0.07 * rs), 6) *
                           ((-1.07 * rs / pow(r + 0.07 * rs, 2)) * (1.12 * pow(rs, 7) / (pow(r, 7) + 0.12 * pow(rs, 7)) - 2) +
                           (-1.12 * pow(rs, 7) * pow(r, 6) / pow(pow(r, 7) + 0.12 * pow(rs, 7), 2)) * (1.07 * rs / (r + 0.07 * rs)));

    QVector<chemkit::Vector3> gradient = distanceGradient(a, b);

    gradient[0] *= de_dr;
    gradient[1] *= de_dr;

    return gradient;
}

// === MmffElectrostaticCalculation ======================================== //
MmffElectrostaticCalculation::MmffElectrostaticCalculation(const MmffAtom *a,
                                                           const MmffAtom *b)
    : MmffCalculation(Electrostatic, 2, 3)
{
    setAtom(0, a);
    setAtom(1, b);
}

bool MmffElectrostaticCalculation::setup(const MmffParameters *parameters)
{
    Q_UNUSED(parameters);

    const MmffAtom *a = atom(0);
    const MmffAtom *b = atom(1);

    chemkit::Float oneFourScaling;

    if(a->isOneFour(b)){
        oneFourScaling = 0.75;
    }
    else{
        oneFourScaling = 1.0;
    }

    setParameter(0, a->charge());
    setParameter(1, b->charge());
    setParameter(2, oneFourScaling);

    return true;
}

chemkit::Float MmffElectrostaticCalculation::energy() const
{
    const chemkit::ForceFieldAtom *a = atom(0);
    const chemkit::ForceFieldAtom *b = atom(1);

    chemkit::Float qa = parameter(0);
    chemkit::Float qb = parameter(1);
    chemkit::Float oneFourScaling = parameter(2);

    chemkit::Float r = distance(a, b);
    chemkit::Float e = 1.0; // dielectric constant
    chemkit::Float d = 0.05; // electrostatic buffering constant

    // equation 13
    return ((332.0716 * qa * qb) / (e * (r + d))) * oneFourScaling;
}

QVector<chemkit::Vector3> MmffElectrostaticCalculation::gradient() const
{
    const chemkit::ForceFieldAtom *a = atom(0);
    const chemkit::ForceFieldAtom *b = atom(1);

    chemkit::Float qa = parameter(0);
    chemkit::Float qb = parameter(1);
    chemkit::Float oneFourScaling = parameter(2);

    chemkit::Float r = distance(a, b);
    chemkit::Float e = 1.0; // dielectric constant
    chemkit::Float d = 0.05; // electrostatic buffering constant

    chemkit::Float de_dr = 332.0716 * qa * qb * oneFourScaling * (-1.0 / (e * pow(r + d, 2)));

    QVector<chemkit::Vector3> gradient = distanceGradient(a, b);

    gradient[0] *= de_dr;
    gradient[1] *= de_dr;

    return gradient;
}
