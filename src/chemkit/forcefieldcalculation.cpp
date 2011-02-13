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

#include "forcefieldcalculation.h"

#include "forcefieldatom.h"

namespace chemkit {

// === ForceFieldCalculationPrivate ======================================== //
class ForceFieldCalculationPrivate
{
    public:
        int type;
        bool setup;
        QVector<Float> parameters;
        QVector<const ForceFieldAtom *> atoms;
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
ForceField* ForceFieldCalculation::forceField()
{
    return const_cast<ForceFieldAtom *>(d->atoms[0])->forceField();
}

/// \overload
const ForceField* ForceFieldCalculation::forceField() const
{
    return d->atoms[0]->forceField();
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
    return d->atoms.value(index, 0);
}

/// Returns the atoms in the calculation.
QVector<const ForceFieldAtom *> ForceFieldCalculation::atoms() const
{
    QVector<const ForceFieldAtom *> atoms(d->atoms.size());

    for(int i = 0; i < d->atoms.size(); i++){
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
    return d->atoms.contains(atom);
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
    return d->parameters.value(index, 0);
}

/// Returns all of the parameters in the calculation.
QVector<Float> ForceFieldCalculation::parameters() const
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
QVector<Vector> ForceFieldCalculation::gradient() const
{
    return numericalGradient();
}

/// Returns the gradient of the energy with respect to the
/// coordinates of each atom for the calculation. This method
/// is used when analytical gradients are not available.
///
/// \see ForceFieldCalculation::gradient()
QVector<Vector> ForceFieldCalculation::numericalGradient() const
{
    QVector<Vector> gradient(atomCount());

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

        gradient[i] = Vector(dx, dy, dz);
    }

    return gradient;
}

// --- Internal Methods ---------------------------------------------------- //
void ForceFieldCalculation::setSetup(bool setup)
{
    d->setup = setup;
}

} // end chemkit namespace
