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

#include <chemkit/cartesiancoordinates.h>

#include "topology.h"
#include "forcefield.h"

namespace chemkit {

// === ForceFieldCalculationPrivate ======================================== //
class ForceFieldCalculationPrivate
{
public:
    ForceField *forceField;
    int type;
    bool setup;
    std::vector<Real> parameters;
    std::vector<size_t> atoms;
};

// === ForceFieldCalculation =============================================== //
/// \class ForceFieldCalculation forcefieldcalculation.h chemkit/forcefieldcalculation.h
/// \ingroup chemkit-md
/// \brief The ForceFieldCalculation class represents an energy
///        calculation in a force field.
///
/// \see ForceField

// --- Construction and Destruction ---------------------------------------- //
ForceFieldCalculation::ForceFieldCalculation(int type,
                                             size_t atomCount,
                                             size_t parameterCount)
    : d(new ForceFieldCalculationPrivate)
{
    d->forceField = 0;
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
    return d->forceField;
}

boost::shared_ptr<Topology> ForceFieldCalculation::topology() const
{
    return d->forceField->topology();
}

// --- Atoms --------------------------------------------------------------- //
/// Sets the atom at \p index to \p atom.
void ForceFieldCalculation::setAtom(size_t index, size_t atom)
{
    d->atoms[index] = atom;
}

/// Returns the atom at index in the calculation.
size_t ForceFieldCalculation::atom(size_t index) const
{
    return d->atoms[index];
}

/// Returns the atoms in the calculation.
std::vector<size_t> ForceFieldCalculation::atoms() const
{
    return d->atoms;
}

/// Returns the number of atoms in the calculation.
size_t ForceFieldCalculation::atomCount() const
{
    return d->atoms.size();
}

std::string ForceFieldCalculation::atomType(size_t index) const
{
    return topology()->type(atom(index));
}

// --- Parameters ---------------------------------------------------------- //
/// Sets the parameter at index to value.
void ForceFieldCalculation::setParameter(int index, Real value)
{
    d->parameters[index] = value;
}

/// Returns the parameter at index.
Real ForceFieldCalculation::parameter(int index) const
{
    return d->parameters[index];
}

/// Returns all of the parameters in the calculation.
std::vector<Real> ForceFieldCalculation::parameters() const
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
Real ForceFieldCalculation::energy(const CartesianCoordinates *coordinates) const
{
    CHEMKIT_UNUSED(coordinates);

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
std::vector<Vector3> ForceFieldCalculation::gradient(const CartesianCoordinates *coordinates) const
{
    return numericalGradient(coordinates);
}

/// Returns the gradient of the energy with respect to the
/// coordinates of each atom for the calculation. This method
/// is used when analytical gradients are not available.
///
/// \see ForceFieldCalculation::gradient()
std::vector<Vector3> ForceFieldCalculation::numericalGradient(const CartesianCoordinates *coordinates) const
{
    std::vector<Vector3> gradient(atomCount());

    CartesianCoordinates writeableCoordinates = *coordinates;

    for(size_t i = 0; i < atomCount(); i++){
        size_t atom = this->atom(i);
        const Point3 &position = coordinates->position(atom);

        // initial energy
        Real eI = energy(&writeableCoordinates);
        Real epsilon = 1.0e-10;

        writeableCoordinates.setPosition(atom, position + Vector3(epsilon, 0, 0));
        Real eF_x = energy(&writeableCoordinates);

        writeableCoordinates.setPosition(atom, position + Vector3(0, epsilon, 0));
        Real eF_y = energy(&writeableCoordinates);

        writeableCoordinates.setPosition(atom, position + Vector3(0, 0, epsilon));
        Real eF_z = energy(&writeableCoordinates);

        // restore initial position
        writeableCoordinates.setPosition(atom, coordinates->position(atom));

        Real dx = (eF_x - eI) / epsilon;
        Real dy = (eF_y - eI) / epsilon;
        Real dz = (eF_z - eI) / epsilon;

        gradient[i] = Vector3(dx, dy, dz);
    }

    return gradient;
}

// --- Internal Methods ---------------------------------------------------- //
void ForceFieldCalculation::setSetup(bool setup)
{
    d->setup = setup;
}

void ForceFieldCalculation::setForceField(ForceField *forceField)
{
    d->forceField = forceField;
}

} // end chemkit namespace
