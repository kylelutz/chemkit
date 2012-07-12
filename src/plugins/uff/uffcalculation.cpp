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

#include "uffcalculation.h"

#include <boost/algorithm/string.hpp>

#include "uffatomtyper.h"
#include "uffforcefield.h"
#include "uffparameters.h"

#include <chemkit/topology.h>
#include <chemkit/constants.h>
#include <chemkit/cartesiancoordinates.h>

// === UffCalculation ====================================================== //
UffCalculation::UffCalculation(int type, int atomCount, int parameterCount)
    : ForceFieldCalculation(type, atomCount, parameterCount)
{
}

// Returns the parameters for the given atom.
const UffAtomParameters* UffCalculation::parameters(const std::string &type) const
{
    const UffForceField *forceField = static_cast<const UffForceField *>(this->forceField());

    return forceField->parameters()->parameters(type);
}

// Returns the bond order of the bond between atom's a and b. If both
// atoms have a resonant type the bond order returned is 1.5.
// Otherwise the integer value of the bond order is returned.
chemkit::Real UffCalculation::bondOrder(size_t a, size_t b) const
{
    const boost::shared_ptr<chemkit::Topology> &topology = this->topology();

    int type = topology->bondedInteractionType(a, b);

    if(type == UffAtomTyper::Resonant){
        return 1.5; // resonant
    }
    else{
        return type;
    }
}

// Returns the length of the bond between two atoms.
chemkit::Real UffCalculation::bondLength(const UffAtomParameters *a, const UffAtomParameters *b, chemkit::Real bondOrder) const
{
    // r_ij = r_i + r_j + r_bo - r_en
    chemkit::Real r_bo = -0.1332 * (a->r + b->r) * log(bondOrder);
    chemkit::Real r_en = ((a->r * b->r) * pow((sqrt(a->X) - sqrt(b->X)), 2)) / (a->X*a->r + b->X*b->r);

    chemkit::Real r_ij = a->r + b->r + r_bo - r_en;

    return r_ij;
}

// === UffBondStrechCalculation ============================================ //
UffBondStrechCalculation::UffBondStrechCalculation(size_t a, size_t b)
    : UffCalculation(BondStrech, 2, 2)
{
    setAtom(0, a);
    setAtom(1, b);
}

bool UffBondStrechCalculation::setup()
{
    const UffAtomParameters *pa = parameters(atomType(0));
    const UffAtomParameters *pb = parameters(atomType(1));

    if(!pa || !pb){
        return false;
    }

    // n = bondorder (1.5 for aromatic, 1.366 for amide)
    chemkit::Real bondorder = bondOrder(atom(0), atom(1));

    chemkit::Real r0 = bondLength(pa, pb, bondorder);

    // parameter(1) = k_ij = 664.12 * (Z*_i * Z*_j) / r_ij^3
    chemkit::Real za = pa->Z;
    chemkit::Real zb = pb->Z;
    chemkit::Real kb = 664.12 * (za * zb) / pow(r0, 3);

    setParameter(0, kb);
    setParameter(1, r0);

    return true;
}

chemkit::Real UffBondStrechCalculation::energy(const chemkit::CartesianCoordinates *coordinates) const
{
    size_t a = atom(0);
    size_t b = atom(1);

    chemkit::Real kb = parameter(0);
    chemkit::Real r0 = parameter(1);
    chemkit::Real r = coordinates->distance(a, b);

    return 0.5 * kb * pow(r - r0, 2);
}

std::vector<chemkit::Vector3> UffBondStrechCalculation::gradient(const chemkit::CartesianCoordinates *coordinates) const
{
    size_t a = atom(0);
    size_t b = atom(1);

    chemkit::Real kb = parameter(0);
    chemkit::Real r0 = parameter(1);
    chemkit::Real r = coordinates->distance(a, b);

    // dE/dr
    chemkit::Real de_dr = kb * (r - r0);

    boost::array<chemkit::Vector3, 2> gradient = coordinates->distanceGradient(a, b);

    gradient[0] *= de_dr;
    gradient[1] *= de_dr;

    return std::vector<chemkit::Vector3>(gradient.begin(), gradient.end());
}

// === UffAngleBendCalculation ============================================= //
UffAngleBendCalculation::UffAngleBendCalculation(size_t a, size_t b, size_t c)
    : UffCalculation(AngleBend, 3, 4)
{
    setAtom(0, a);
    setAtom(1, b);
    setAtom(2, c);
}

bool UffAngleBendCalculation::setup()
{
    const UffAtomParameters *pa = parameters(atomType(0));
    const UffAtomParameters *pb = parameters(atomType(1));
    const UffAtomParameters *pc = parameters(atomType(2));

    if(!pa || !pb || !pc){
        return false;
    }

    chemkit::Real theta0 = pb->theta * chemkit::constants::DegreesToRadians;

    chemkit::Real bo_ij = bondOrder(atom(0), atom(1));
    chemkit::Real bo_jk = bondOrder(atom(1), atom(2));

    chemkit::Real r_ab = bondLength(pa, pb, bo_ij);
    chemkit::Real r_bc = bondLength(pb, pc, bo_jk);
    chemkit::Real r_ac = sqrt(pow(r_ab, 2)  + pow(r_bc, 2) - (2.0 * r_ab * r_bc * cos(theta0)));

    chemkit::Real beta = 664.12 / (r_ab * r_bc);

    chemkit::Real z_a = pa->Z;
    chemkit::Real z_c = pc->Z;

    // equation 13
    chemkit::Real ka = beta * ((z_a * z_c) / pow(r_ac, 5)) * r_ab * r_bc * (3.0 * r_ab * r_bc * (1.0 - pow(cos(theta0), 2.0)) - (pow(r_ac, 2.0) * cos(theta0)));

    setParameter(0, ka);

    chemkit::Real sinTheta0 = sin(theta0);

    // clamp sin(theta0) to 1e-3 because for some atoms theta0 == pi which
    // would lead to a division by zero error when calculating c2 below
    if(std::abs(sinTheta0) < 1e-3){
        sinTheta0 = 1e-3;
    }

    chemkit::Real c2 = 1 / (4 * pow(sinTheta0, 2));
    chemkit::Real c1 = -4 * c2 * cos(theta0);
    chemkit::Real c0 = c2 * (2 * pow(cos(theta0), 2) + 1);

    setParameter(1, c0);
    setParameter(2, c1);
    setParameter(3, c2);

    return true;
}

chemkit::Real UffAngleBendCalculation::energy(const chemkit::CartesianCoordinates *coordinates) const
{
    size_t a = atom(0);
    size_t b = atom(1);
    size_t c = atom(2);

    chemkit::Real ka = parameter(0);
    chemkit::Real c0 = parameter(1);
    chemkit::Real c1 = parameter(2);
    chemkit::Real c2 = parameter(3);

    chemkit::Real theta = coordinates->angleRadians(a, b, c);

    return ka * (c0 + (c1 * cos(theta)) + (c2 * cos(2*theta)));
}

std::vector<chemkit::Vector3> UffAngleBendCalculation::gradient(const chemkit::CartesianCoordinates *coordinates) const
{
    size_t a = atom(0);
    size_t b = atom(1);
    size_t c = atom(2);

    chemkit::Real ka = parameter(0);
    chemkit::Real c1 = parameter(2);
    chemkit::Real c2 = parameter(3);

    chemkit::Real theta = coordinates->angleRadians(a, b, c);

    // dE/dtheta
    chemkit::Real de_dtheta = -ka * (c1 * sin(theta) + 2 * c2 * sin(2 * theta));

    boost::array<chemkit::Vector3, 3> gradient = coordinates->angleGradientRadians(a, b, c);

    gradient[0] *= de_dtheta;
    gradient[1] *= de_dtheta;
    gradient[2] *= de_dtheta;

    return std::vector<chemkit::Vector3>(gradient.begin(), gradient.end());
}

// === UffTorsionCalculation =============================================== //
UffTorsionCalculation::UffTorsionCalculation(size_t a, size_t b, size_t c, size_t d)
    : UffCalculation(Torsion, 4, 3)
{
    setAtom(0, a);
    setAtom(1, b);
    setAtom(2, c);
    setAtom(3, d);
}

bool UffTorsionCalculation::setup()
{
    UffForceField *forceField = static_cast<UffForceField *>(this->forceField());
    const boost::shared_ptr<chemkit::Topology> &topology = this->topology();

    size_t b = atom(1);
    size_t c = atom(2);

    std::string typeB = topology->type(b);
    std::string typeC = topology->type(c);

    if(typeB.length() < 3 || typeC.length() < 3){
        return false;
    }

    const UffAtomParameters *pb = parameters(typeB);
    const UffAtomParameters *pc = parameters(typeC);

    chemkit::Real V = 0;
    chemkit::Real n = 0;
    chemkit::Real phi0 = 0;

    // sp3-sp3
    if(typeB[2] == '3' && typeC[2] == '3'){

        // exception for two group six atoms
        if(forceField->isGroupSix(atom(1)) && forceField->isGroupSix(atom(2))){
            if(boost::starts_with(typeB, "O_") && boost::starts_with(typeC, "O_")){
                V = 2; // sqrt(2*2)
            }
            else if(boost::starts_with(typeB, "O_") || boost::starts_with(typeC, "O_")){
                V = sqrt(2 * 6.8);
            }
            else{
                V = sqrt(6.8 * 6.8);
            }

            n = 2;
            phi0 = 90;
        }

        // general case
        else{
            // equation 16
            V = sqrt(pb->V * pc->V);

            n = 3;
            phi0 = 180 * chemkit::constants::DegreesToRadians;
        }
    }
    // sp2-sp2
    else if((typeB[2] == '2' || typeB[2] == 'R') && (typeC[2] == '2' || typeC[2] == 'R')){
        chemkit::Real bondorder = bondOrder(b, c);

        // equation 17
        V = 5 * sqrt(pb->U * pc->U) * (1 + 4.18 * log(bondorder));

        n = 2;
        phi0 = 180 * chemkit::constants::DegreesToRadians;
    }
    // group 6 sp3 - any sp2 or R
    else if((forceField->isGroupSix(b) && (typeC[2] == '2' || typeC[2] == 'R')) ||
            (forceField->isGroupSix(c) && (typeB[2] == '2' || typeB[2] == 'R'))){
        chemkit::Real bondorder = bondOrder(b, c);

        // equation 17
        V = 5 * sqrt(pb->U * pc->U) * (1 + 4.18 * log(bondorder));

        n = 2;
        phi0 = 90 * chemkit::constants::DegreesToRadians;
    }
    // sp3-sp2
    else if((typeB[2] == '3' && (typeC[2] == '2' || typeC[2] == 'R')) ||
            (typeC[2] == '3' && (typeB[2] == '2' || typeB[2] == 'R'))){
        V = 1;
        n = 6;
        phi0 = 0;
    }
    else{
        return false;
    }

    setParameter(0, V);
    setParameter(1, n);
    setParameter(2, phi0);

    return true;
}

chemkit::Real UffTorsionCalculation::energy(const chemkit::CartesianCoordinates *coordinates) const
{
    size_t a = atom(0);
    size_t b = atom(1);
    size_t c = atom(2);
    size_t d = atom(3);

    chemkit::Real V = parameter(0);
    chemkit::Real n = parameter(1);
    chemkit::Real phi0 = parameter(2);

    chemkit::Real phi = coordinates->torsionAngleRadians(a, b, c, d);

    return 0.5 * V * (1 - cos(n * phi0) * cos(n * phi));
}

std::vector<chemkit::Vector3> UffTorsionCalculation::gradient(const chemkit::CartesianCoordinates *coordinates) const
{
    size_t a = atom(0);
    size_t b = atom(1);
    size_t c = atom(2);
    size_t d = atom(3);

    chemkit::Real V = parameter(0);
    chemkit::Real n = parameter(1);
    chemkit::Real phi0 = parameter(2);

    chemkit::Real phi = coordinates->torsionAngleRadians(a, b, c, d);

    // dE/dphi
    chemkit::Real de_dphi = 0.5 * V * n * cos(n * phi0) * sin(n * phi);

    boost::array<chemkit::Vector3, 4> gradient = coordinates->torsionAngleGradientRadians(a, b, c, d);

    gradient[0] *= de_dphi;
    gradient[1] *= de_dphi;
    gradient[2] *= de_dphi;
    gradient[3] *= de_dphi;

    return std::vector<chemkit::Vector3>(gradient.begin(), gradient.end());
}

// === UffInversionCalculation ============================================= //
UffInversionCalculation::UffInversionCalculation(size_t a, size_t b, size_t c, size_t d)
    : UffCalculation(Inversion, 4, 4)
{
    setAtom(0, a);
    setAtom(1, b);
    setAtom(2, c);
    setAtom(3, d);
}

bool UffInversionCalculation::setup()
{
    const boost::shared_ptr<chemkit::Topology> &topology = this->topology();

    // b is the center atom
    size_t a = atom(0);
    size_t b = atom(1);
    size_t c = atom(2);
    size_t d = atom(3);

    std::string typeA = topology->type(a);
    std::string typeB = topology->type(b);
    std::string typeC = topology->type(c);
    std::string typeD = topology->type(d);


    chemkit::Real k = 0;
    chemkit::Real c0 = 0;
    chemkit::Real c1 = 0;
    chemkit::Real c2 = 0;

    // sp2 carbon
    if(typeB == "C_2" || typeB == "C_R"){
        if(typeA == "O_2" || typeC == "O_2" || typeD == "O_2"){
            k = 50;
        }
        else{
            k = 6;
        }

        c0 = 1;
        c1 = -1;
        c2 = 0;
    }

    // divide by 3
    k /= 3;

    setParameter(0, k);
    setParameter(1, c0);
    setParameter(2, c1);
    setParameter(3, c2);

    return true;
}

chemkit::Real UffInversionCalculation::energy(const chemkit::CartesianCoordinates *coordinates) const
{
    size_t a = atom(0);
    size_t b = atom(1);
    size_t c = atom(2);
    size_t d = atom(3);

    chemkit::Real k = parameter(0);
    chemkit::Real c0 = parameter(1);
    chemkit::Real c1 = parameter(2);
    chemkit::Real c2 = parameter(3);

    chemkit::Real w = coordinates->wilsonAngleRadians(a, b, c, d);
    chemkit::Real y = w + (chemkit::constants::Pi / 2.0);

    return k * (c0 + c1 * sin(y) + c2 * cos(2 * y));
}

std::vector<chemkit::Vector3> UffInversionCalculation::gradient(const chemkit::CartesianCoordinates *coordinates) const
{
    size_t a = atom(0);
    size_t b = atom(1);
    size_t c = atom(2);
    size_t d = atom(3);

    chemkit::Real k = parameter(0);
    chemkit::Real c1 = parameter(2);
    chemkit::Real c2 = parameter(3);

    chemkit::Real w = coordinates->wilsonAngleRadians(a, b, c, d);
    chemkit::Real y = w + (chemkit::constants::Pi / 2.0);

    // dE/dw
    chemkit::Real de_dw = k * (c1 * cos(y) - 2 * c2 * sin(2 * y));

    boost::array<chemkit::Vector3, 4> gradient = coordinates->wilsonAngleGradientRadians(a, b, c, d);

    gradient[0] *= de_dw;
    gradient[1] *= de_dw;
    gradient[2] *= de_dw;
    gradient[3] *= de_dw;

    return std::vector<chemkit::Vector3>(gradient.begin(), gradient.end());
}

// === UffVanDerWaalsCalculation =========================================== //
UffVanDerWaalsCalculation::UffVanDerWaalsCalculation(size_t a, size_t b)
    : UffCalculation(VanDerWaals, 2, 2)
{
    setAtom(0, a);
    setAtom(1, b);
}

bool UffVanDerWaalsCalculation::setup()
{
    const UffAtomParameters *pa = parameters(atomType(0));
    const UffAtomParameters *pb = parameters(atomType(1));
    if(!pa || !pb){
        return false;
    }

    // equation 22
    chemkit::Real d = sqrt(pa->D * pb->D);

    // equation 21b
    chemkit::Real x = sqrt(pa->x * pb->x);

    setParameter(0, d);
    setParameter(1, x);

    return true;
}

chemkit::Real UffVanDerWaalsCalculation::energy(const chemkit::CartesianCoordinates *coordinates) const
{
    size_t a = atom(0);
    size_t b = atom(1);

    chemkit::Real d = parameter(0);
    chemkit::Real x = parameter(1);
    chemkit::Real r = coordinates->distance(a, b);

    return d * (-2 * pow(x/r, 6) + pow(x/r, 12));
}

std::vector<chemkit::Vector3> UffVanDerWaalsCalculation::gradient(const chemkit::CartesianCoordinates *coordinates) const
{
    size_t a = atom(0);
    size_t b = atom(1);

    chemkit::Real d = parameter(0);
    chemkit::Real x = parameter(1);
    chemkit::Real r = coordinates->distance(a, b);

    // dE/dr
    chemkit::Real de_dr = -12 * d * x / pow(r, 2) * (pow(x/r, 11) - pow(x/r, 5));

    boost::array<chemkit::Vector3, 2> gradient = coordinates->distanceGradient(a, b);

    gradient[0] *= de_dr;
    gradient[1] *= de_dr;

    return std::vector<chemkit::Vector3>(gradient.begin(), gradient.end());
}

// === UffElectrostaticCalculation ========================================= //
UffElectrostaticCalculation::UffElectrostaticCalculation(size_t a, size_t b)
    : UffCalculation(Electrostatic, 2, 2)
{
    setAtom(0, a);
    setAtom(1, b);
}

bool UffElectrostaticCalculation::setup()
{
    return false;
}

chemkit::Real UffElectrostaticCalculation::energy(const chemkit::CartesianCoordinates *coordinates) const
{
    size_t a = atom(0);
    size_t b = atom(1);

    chemkit::Real qa = parameter(0);
    chemkit::Real qb = parameter(1);

    chemkit::Real e = 1;
    chemkit::Real r = coordinates->distance(a, b);

    return 332.037 * (qa * qb) / (e * r);
}
