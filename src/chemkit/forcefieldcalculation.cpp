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

#include "forcefieldcalculation.h"

#include <algorithm>

#include "forcefieldatom.h"

namespace chemkit {

// === ForceFieldCalculationPrivate ======================================== //
class ForceFieldCalculationPrivate
{
    public:
        int type;
        bool setup;
        std::vector<Float> parameters;
        std::vector<const ForceFieldAtom *> atoms;
};

// === ForceFieldCalculation =============================================== //
/// \class ForceFieldCalculation forcefieldcalculation.h chemkit/forcefieldcalculation.h
/// \ingroup chemkit
/// \brief The ForceFieldCalculation class represents an energy
///        calculation in a force field.
///
/// \see ForceField

// --- Construction and Destruction ---------------------------------------- //
ForceFieldCalculation::ForceFieldCalculation(int type, int atomCount, int parameterCount)
    : d(new ForceFieldCalculationPrivate)
{
    d->type = type;
    d->setup = false;
    d->atoms.resize(atomCount);
    d->parameters.resize(parameterCount);
}

ForceFieldCalculation::~ForceFieldCalculation()
{
    delete d;
}

// --- Properties ---------------------------------------------------------- //
/// Returns the type of the calculation.
int ForceFieldCalculation::type() const
{
    return d->type;
}

/// Returns \c true if the calculation is setup.
bool ForceFieldCalculation::isSetup() const
{
    return d->setup;
}

/// Returns the force field the calculation is a part of.
ForceField* ForceFieldCalculation::forceField() const
{
    return const_cast<ForceFieldAtom *>(d->atoms[0])->forceField();
}

// --- Atoms --------------------------------------------------------------- //
/// Sets the atom at \p index to \p atom.
void ForceFieldCalculation::setAtom(int index, const ForceFieldAtom *atom)
{
    d->atoms[index] = atom;
}

/// Returns the atom at index in the calculation.
const ForceFieldAtom* ForceFieldCalculation::atom(int index) const
{
    return d->atoms[index];
}

/// Returns the atoms in the calculation.
std::vector<const ForceFieldAtom *> ForceFieldCalculation::atoms() const
{
    std::vector<const ForceFieldAtom *> atoms(d->atoms.size());

    for(unsigned int i = 0; i < d->atoms.size(); i++){
        atoms[i] = d->atoms[i];
    }

    return atoms;
}

/// Returns the number of atoms in the calculation.
int ForceFieldCalculation::atomCount() const
{
    return d->atoms.size();
}

/// Returns \c true if the calculation contains the atom.
bool ForceFieldCalculation::contains(const ForceFieldAtom *atom) const
{
    return std::find(d->atoms.begin(), d->atoms.end(), atom) != d->atoms.end();
}

// --- Parameters ---------------------------------------------------------- //
/// Sets the parameter at index to value.
void ForceFieldCalculation::setParameter(int index, Float value)
{
    d->parameters[index] = value;
}

/// Returns the parameter at index.
Float ForceFieldCalculation::parameter(int index) const
{
    return d->parameters[index];
}

/// Returns all of the parameters in the calculation.
std::vector<Float> ForceFieldCalculation::parameters() const
{
    return d->parameters;
}

/// Returns the number of parameters in the calculation.
int ForceFieldCalculation::parameterCount() const
{
    return d->parameters.size();
}

// --- Calculations -------------------------------------------------------- //
/// Returns the energy of the calculation. Energy is in kcal/mol.
Float ForceFieldCalculation::energy() const
{
    return 0;
}

/// Returns the gradient of the energy with respect to the
/// coordinates of each atom for the calculation. This method will
/// calculate the derivative analytically if possible, otherwise it
/// will be calculated numerically via numericalGradient().
///
/** \f[ \nabla E = \left[
///                \begin{array}{ccc}
///                    \frac{\partial E}{\partial x_{0}} &
///                    \frac{\partial E}{\partial y_{0}} &
///                    \frac{\partial E}{\partial z_{0}} \\
///                    \frac{\partial E}{\partial x_{1}} &
///                    \frac{\partial E}{\partial y_{1}} &
///                    \frac{\partial E}{\partial z_{1}} \\
///                    \vdots & \vdots & \vdots \\
///                    \frac{\partial E}{\partial x_{n}} &
///                    \frac{\partial E}{\partial y_{n}} &
///                    \frac{\partial E}{\partial z_{n}}
///                \end{array}
///                \right]
/// \f]
**/
std::vector<Vector3> ForceFieldCalculation::gradient() const
{
    return numericalGradient();
}

/// Returns the gradient of the energy with respect to the
/// coordinates of each atom for the calculation. This method
/// is used when analytical gradients are not available.
///
/// \see ForceFieldCalculation::gradient()
std::vector<Vector3> ForceFieldCalculation::numericalGradient() const
{
    std::vector<Vector3> gradient(atomCount());

    for(int i = 0; i < atomCount(); i++){
        ForceFieldAtom *atom = const_cast<ForceFieldAtom *>(d->atoms[i]);

        // initial energy
        Float eI = energy();
        Float epsilon = 1.0e-10;

        atom->moveBy(epsilon, 0, 0);
        Float eF_x = energy();

        atom->moveBy(-epsilon, epsilon, 0);
        Float eF_y = energy();

        atom->moveBy(0, -epsilon, epsilon);
        Float eF_z = energy();

        // restore initial position
        atom->moveBy(0, 0, -epsilon);

        Float dx = (eF_x - eI) / epsilon;
        Float dy = (eF_y - eI) / epsilon;
        Float dz = (eF_z - eI) / epsilon;

        gradient[i] = Vector3(dx, dy, dz);
    }

    return gradient;
}

// --- Internal Methods ---------------------------------------------------- //
void ForceFieldCalculation::setSetup(bool setup)
{
    d->setup = setup;
}

} // end chemkit namespace
