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
#include "uffforcefield.h"
#include "uffparameters.h"

#include <chemkit/atom.h>
#include <chemkit/bond.h>
#include <chemkit/molecule.h>
#include <chemkit/constants.h>

// === UffCalculation ====================================================== //
UffCalculation::UffCalculation(int type, int atomCount, int parameterCount)
    : ForceFieldCalculation(type, atomCount, parameterCount)
{
}

// Returns the parameters for the given atom.
const UffAtomParameters* UffCalculation::parameters(const chemkit::ForceFieldAtom *atom) const
{
    const UffForceField *forceField = static_cast<const UffForceField *>(this->forceField());
    const UffParameters *parameters = forceField->parameters();

    return parameters->parameters(atom);
}

// Returns the bond order of the bond between atom's a and b. If both
// atoms have a resonant type the bond order returned is 1.5.
// Otherwise the integer value of the bond order is returned.
chemkit::Real UffCalculation::bondOrder(const chemkit::ForceFieldAtom *a, const chemkit::ForceFieldAtom *b) const
{
    const chemkit::Bond *bond = a->atom()->bondTo(b->atom());

    if((a->type().length() > 2 && a->type()[2] == 'R') && (b->type().length() > 2 && b->type()[2] == 'R')){
        return 1.5; // resonant
    }
    else{
        return bond->order();
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
UffBondStrechCalculation::UffBondStrechCalculation(const chemkit::ForceFieldAtom *a, const chemkit::ForceFieldAtom *b)
    : UffCalculation(BondStrech, 2, 2)
{
    setAtom(0, a);
    setAtom(1, b);
}

bool UffBondStrechCalculation::setup()
{
    const UffAtomParameters *pa = parameters(atom(0));
    const UffAtomParameters *pb = parameters(atom(1));

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

chemkit::Real UffBondStrechCalculation::energy() const
{
    const chemkit::ForceFieldAtom *a = atom(0);
    const chemkit::ForceFieldAtom *b = atom(1);

    chemkit::Real kb = parameter(0);
    chemkit::Real r0 = parameter(1);
    chemkit::Real r = distance(a, b);

    return 0.5 * kb * pow(r - r0, 2);
}

std::vector<chemkit::Vector3> UffBondStrechCalculation::gradient() const
{
    const chemkit::ForceFieldAtom *a = atom(0);
    const chemkit::ForceFieldAtom *b = atom(1);

    chemkit::Real kb = parameter(0);
    chemkit::Real r0 = parameter(1);
    chemkit::Real r = distance(a, b);

    // dE/dr
    chemkit::Real de_dr = kb * (r - r0);

    boost::array<chemkit::Vector3, 2> gradient = distanceGradient(a, b);

    gradient[0] *= de_dr;
    gradient[1] *= de_dr;

    return std::vector<chemkit::Vector3>(gradient.begin(), gradient.end());
}

// === UffAngleBendCalculation ============================================= //
UffAngleBendCalculation::UffAngleBendCalculation(const chemkit::ForceFieldAtom *a,
                                                 const chemkit::ForceFieldAtom *b,
                                                 const chemkit::ForceFieldAtom *c)
    : UffCalculation(AngleBend, 3, 4)
{
    setAtom(0, a);
    setAtom(1, b);
    setAtom(2, c);
}

bool UffAngleBendCalculation::setup()
{
    const UffAtomParameters *pa = parameters(atom(0));
    const UffAtomParameters *pb = parameters(atom(1));
    const UffAtomParameters *pc = parameters(atom(2));

    if(!pa || !pb || !pc){
        return false;
    }

    chemkit::Real theta0 = pb->theta * chemkit::constants::DegreesToRadians;

    const chemkit::Bond *bond_ab = atom(0)->atom()->bondTo(atom(1)->atom());
    const chemkit::Bond *bond_bc = atom(1)->atom()->bondTo(atom(2)->atom());

    chemkit::Real bo_ij = bond_ab->order();
    chemkit::Real bo_jk = bond_bc->order();

    chemkit::Real r_ab = bondLength(pa, pb, bo_ij);
    chemkit::Real r_bc = bondLength(pb, pc, bo_jk);
    chemkit::Real r_ac = sqrt(pow(r_ab, 2)  + pow(r_bc, 2) - (2.0 * r_ab * r_bc * cos(theta0)));

    chemkit::Real beta = 664.12 / (r_ab * r_bc);

    chemkit::Real z_a = pa->Z;
    chemkit::Real z_c = pc->Z;

    // equation 13
    chemkit::Real ka = beta * ((z_a * z_c) / pow(r_ac, 5)) * r_ab * r_bc * (3.0 * r_ab * r_bc * (1.0 - pow(cos(theta0), 2.0)) - (pow(r_ac, 2.0) * cos(theta0)));

    setParameter(0, ka);

    chemkit::Real c2 = 1 / (4 * pow(sin(theta0), 2));
    chemkit::Real c1 = -4 * c2 * cos(theta0);
    chemkit::Real c0 = c2 * (2 * pow(cos(theta0), 2) + 1);

    setParameter(1, c0);
    setParameter(2, c1);
    setParameter(3, c2);

    return true;
}

chemkit::Real UffAngleBendCalculation::energy() const
{
    const chemkit::ForceFieldAtom *a = atom(0);
    const chemkit::ForceFieldAtom *b = atom(1);
    const chemkit::ForceFieldAtom *c = atom(2);

    chemkit::Real ka = parameter(0);
    chemkit::Real c0 = parameter(1);
    chemkit::Real c1 = parameter(2);
    chemkit::Real c2 = parameter(3);

    chemkit::Real theta = bondAngleRadians(a, b, c);

    return ka * (c0 + (c1 * cos(theta)) + (c2 * cos(2*theta)));
}

std::vector<chemkit::Vector3> UffAngleBendCalculation::gradient() const
{
    const chemkit::ForceFieldAtom *a = atom(0);
    const chemkit::ForceFieldAtom *b = atom(1);
    const chemkit::ForceFieldAtom *c = atom(2);

    chemkit::Real ka = parameter(0);
    chemkit::Real c1 = parameter(2);
    chemkit::Real c2 = parameter(3);

    chemkit::Real theta = bondAngleRadians(a, b, c);

    // dE/dtheta
    chemkit::Real de_dtheta = -ka * (c1 * sin(theta) + 2 * c2 * sin(2 * theta));

    boost::array<chemkit::Vector3, 3> gradient = bondAngleGradientRadians(a, b, c);

    gradient[0] *= de_dtheta;
    gradient[1] *= de_dtheta;
    gradient[2] *= de_dtheta;

    return std::vector<chemkit::Vector3>(gradient.begin(), gradient.end());
}

// === UffTorsionCalculation =============================================== //
UffTorsionCalculation::UffTorsionCalculation(const chemkit::ForceFieldAtom *a,
                                             const chemkit::ForceFieldAtom *b,
                                             const chemkit::ForceFieldAtom *c,
                                             const chemkit::ForceFieldAtom *d)
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

    const chemkit::ForceFieldAtom *b = atom(1);
    const chemkit::ForceFieldAtom *c = atom(2);

    if(b->type().length() < 3 || c->type().length() < 3){
        return false;
    }

    const UffAtomParameters *pb = parameters(b);
    const UffAtomParameters *pc = parameters(c);

    chemkit::Real V = 0;
    chemkit::Real n = 0;
    chemkit::Real phi0 = 0;

    // sp3-sp3
    if(b->type()[2] == '3' && c->type()[2] == '3'){

        // exception for two group six atoms
        if(forceField->isGroupSix(b) && forceField->isGroupSix(c)){
            if(b->atom()->is(chemkit::Atom::Oxygen) && c->atom()->is(chemkit::Atom::Oxygen)){
                V = 2; // sqrt(2*2)
            }
            else if(b->atom()->is(chemkit::Atom::Oxygen) || c->atom()->is(chemkit::Atom::Oxygen)){
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
    else if((b->type()[2] == '2' || b->type()[2] == 'R') && (c->type()[2] == '2' || c->type()[2] == 'R')){
        chemkit::Real bondorder = bondOrder(b, c);

        // equation 17
        V = 5 * sqrt(pb->U * pc->U) * (1 + 4.18 * log(bondorder));

        n = 2;
        phi0 = 180 * chemkit::constants::DegreesToRadians;
    }
    // group 6 sp3 - any sp2 or R
    else if((forceField->isGroupSix(b) && (c->type()[2] == '2' || c->type()[2] == 'R')) ||
            (forceField->isGroupSix(c) && (b->type()[2] == '2' || b->type()[2] == 'R'))){
        chemkit::Real bondorder = bondOrder(b, c);

        // equation 17
        V = 5 * sqrt(pb->U * pc->U) * (1 + 4.18 * log(bondorder));

        n = 2;
        phi0 = 90 * chemkit::constants::DegreesToRadians;
    }
    // sp3-sp2
    else if((b->type()[2] == '3' && (c->type()[2] == '2' || c->type()[2] == 'R')) ||
            (c->type()[2] == '3' && (b->type()[2] == '2' || b->type()[2] == 'R'))){
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

chemkit::Real UffTorsionCalculation::energy() const
{
    const chemkit::ForceFieldAtom *a = atom(0);
    const chemkit::ForceFieldAtom *b = atom(1);
    const chemkit::ForceFieldAtom *c = atom(2);
    const chemkit::ForceFieldAtom *d = atom(3);

    chemkit::Real V = parameter(0);
    chemkit::Real n = parameter(1);
    chemkit::Real phi0 = parameter(2);

    chemkit::Real phi = torsionAngleRadians(a, b, c, d);

    return 0.5 * V * (1 - cos(n * phi0) * cos(n * phi));
}

std::vector<chemkit::Vector3> UffTorsionCalculation::gradient() const
{
    const chemkit::ForceFieldAtom *a = atom(0);
    const chemkit::ForceFieldAtom *b = atom(1);
    const chemkit::ForceFieldAtom *c = atom(2);
    const chemkit::ForceFieldAtom *d = atom(3);

    chemkit::Real V = parameter(0);
    chemkit::Real n = parameter(1);
    chemkit::Real phi0 = parameter(2);

    chemkit::Real phi = torsionAngleRadians(a, b, c, d);

    // dE/dphi
    chemkit::Real de_dphi = 0.5 * V * n * cos(n * phi0) * sin(n * phi);

    boost::array<chemkit::Vector3, 4> gradient = torsionAngleGradientRadians(a, b, c, d);

    gradient[0] *= de_dphi;
    gradient[1] *= de_dphi;
    gradient[2] *= de_dphi;
    gradient[3] *= de_dphi;

    return std::vector<chemkit::Vector3>(gradient.begin(), gradient.end());
}

// === UffInversionCalculation ============================================= //
UffInversionCalculation::UffInversionCalculation(const chemkit::ForceFieldAtom *a,
                                                 const chemkit::ForceFieldAtom *b,
                                                 const chemkit::ForceFieldAtom *c,
                                                 const chemkit::ForceFieldAtom *d)
    : UffCalculation(Inversion, 4, 4)
{
    setAtom(0, a);
    setAtom(1, b);
    setAtom(2, c);
    setAtom(3, d);
}

bool UffInversionCalculation::setup()
{
    // b is the center atom
    const chemkit::ForceFieldAtom *a = atom(0);
    const chemkit::ForceFieldAtom *b = atom(1);
    const chemkit::ForceFieldAtom *c = atom(2);
    const chemkit::ForceFieldAtom *d = atom(3);

    chemkit::Real k = 0;
    chemkit::Real c0 = 0;
    chemkit::Real c1 = 0;
    chemkit::Real c2 = 0;

    // sp2 carbon
    if(b->type() == "C_2" || b->type() == "C_R"){
        if(a->type() == "O_2" || c->type() == "O_2" || d->type() == "O_2"){
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

chemkit::Real UffInversionCalculation::energy() const
{
    const chemkit::ForceFieldAtom *a = atom(0);
    const chemkit::ForceFieldAtom *b = atom(1);
    const chemkit::ForceFieldAtom *c = atom(2);
    const chemkit::ForceFieldAtom *d = atom(3);

    chemkit::Real k = parameter(0);
    chemkit::Real c0 = parameter(1);
    chemkit::Real c1 = parameter(2);
    chemkit::Real c2 = parameter(3);

    chemkit::Real w = wilsonAngleRadians(a, b, c, d);
    chemkit::Real y = w + (chemkit::constants::Pi / 2.0);

    return k * (c0 + c1 * sin(y) + c2 * cos(2 * y));
}

std::vector<chemkit::Vector3> UffInversionCalculation::gradient() const
{
    const chemkit::ForceFieldAtom *a = atom(0);
    const chemkit::ForceFieldAtom *b = atom(1);
    const chemkit::ForceFieldAtom *c = atom(2);
    const chemkit::ForceFieldAtom *d = atom(3);

    chemkit::Real k = parameter(0);
    chemkit::Real c1 = parameter(2);
    chemkit::Real c2 = parameter(3);

    chemkit::Real w = wilsonAngleRadians(a, b, c, d);
    chemkit::Real y = w + (chemkit::constants::Pi / 2.0);

    // dE/dw
    chemkit::Real de_dw = k * (c1 * cos(y) - 2 * c2 * sin(2 * y));

    boost::array<chemkit::Vector3, 4> gradient = wilsonAngleGradientRadians(a, b, c, d);

    gradient[0] *= de_dw;
    gradient[1] *= de_dw;
    gradient[2] *= de_dw;
    gradient[3] *= de_dw;

    return std::vector<chemkit::Vector3>(gradient.begin(), gradient.end());
}

// === UffVanDerWaalsCalculation =========================================== //
UffVanDerWaalsCalculation::UffVanDerWaalsCalculation(const chemkit::ForceFieldAtom *a,
                                                     const chemkit::ForceFieldAtom *b)
    : UffCalculation(VanDerWaals, 2, 2)
{
    setAtom(0, a);
    setAtom(1, b);
}

bool UffVanDerWaalsCalculation::setup()
{
    const chemkit::ForceFieldAtom *a = atom(0);
    const chemkit::ForceFieldAtom *b = atom(1);

    const UffAtomParameters *pa = parameters(a);
    const UffAtomParameters *pb = parameters(b);

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

chemkit::Real UffVanDerWaalsCalculation::energy() const
{
    const chemkit::ForceFieldAtom *a = atom(0);
    const chemkit::ForceFieldAtom *b = atom(1);

    chemkit::Real d = parameter(0);
    chemkit::Real x = parameter(1);
    chemkit::Real r = distance(a, b);

    return d * (-2 * pow(x/r, 6) + pow(x/r, 12));
}

std::vector<chemkit::Vector3> UffVanDerWaalsCalculation::gradient() const
{
    const chemkit::ForceFieldAtom *a = atom(0);
    const chemkit::ForceFieldAtom *b = atom(1);

    chemkit::Real d = parameter(0);
    chemkit::Real x = parameter(1);
    chemkit::Real r = distance(a, b);

    // dE/dr
    chemkit::Real de_dr = -12 * d * x / pow(r, 2) * (pow(x/r, 11) - pow(x/r, 5));

    boost::array<chemkit::Vector3, 2> gradient = distanceGradient(a, b);

    gradient[0] *= de_dr;
    gradient[1] *= de_dr;

    return std::vector<chemkit::Vector3>(gradient.begin(), gradient.end());
}

// === UffElectrostaticCalculation ========================================= //
UffElectrostaticCalculation::UffElectrostaticCalculation(const chemkit::ForceFieldAtom *a,
                                                         const chemkit::ForceFieldAtom *b)
    : UffCalculation(Electrostatic, 2, 2)
{
    setAtom(0, a);
    setAtom(1, b);
}

bool UffElectrostaticCalculation::setup()
{
    return false;
}

chemkit::Real UffElectrostaticCalculation::energy() const
{
    const chemkit::ForceFieldAtom *a = atom(0);
    const chemkit::ForceFieldAtom *b = atom(1);

    chemkit::Real qa = parameter(0);
    chemkit::Real qb = parameter(1);

    chemkit::Real e = 1;
    chemkit::Real r = distance(a, b);

    return 332.037 * (qa * qb) / (e * r);
}
