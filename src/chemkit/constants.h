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

#ifndef CHEMKIT_CONSTANTS_H
#define CHEMKIT_CONSTANTS_H

#include "chemkit.h"

namespace chemkit {

/// \ingroup chemkit
/// \brief The %chemkit::%constants namespace contains values for
///        mathematical and physical constants.
namespace constants {

/// Pi (\f$\pi\f$). The ratio of the circumference of a circle
/// to its diameter.
const Float Pi = 3.14159265358979323846264338327;

/// Conversion factor from degrees to radians.
///
/// Equal to \f$\frac{\pi}{180}\f$.
const Float DegreesToRadians = Pi / 180.0;

/// Conversion factor from radians to degrees.
///
/// Equal to \f$\frac{180}{\pi}\f$.
const Float RadiansToDegrees = 180.0 / Pi;

/// Conversion factor from calories to joules.
const Float CaloriesToJoules = 4.184;

/// Conversion factor from joules to calories.
const Float JoulesToCalories = 1.0 / CaloriesToJoules;

/// Avogadro constant (\f$N_{A}\f$).
const Float AvogadroConstant = 6.02214179e23;

/// Boltzmann constant (\f$k_{B}\f$)
const Float BoltzmannConstant = 1.3806503e-23;

/// Gas constant.
///
/// Equal to \f$N_{A} k_{B}\f$.
const Float GasConstant = AvogadroConstant * BoltzmannConstant;

/// Planck constant (\f$h\f$).
const Float PlanckConstant = 6.62606896e-34;

/// Reduced Planck constant (\f$\hbar\f$).
const Float ReducedPlanckConstant = PlanckConstant / (2.0 * Pi);

/// Mass of a proton (\f$m_{p}\f$).
const Float ProtonMass = 1.672621637e-27;

/// Mass of a electron (\f$m_{e}\f$).
const Float ElectronMass = 9.10938215e-31;

/// Elementary charge (\f$e\f$).
const Float ElementaryCharge = 1.602176487e-19;

/// Charge of a proton (\f$+e\f$).
const Float ProtonCharge = ElementaryCharge;

/// Charge of a electron (\f$-e\f$).
const Float ElectronCharge = -ElementaryCharge;

/// Faraday constant (\f$F\f$).
///
/// Equal to \f$e N_{A}\f$.
const Float FaradayConstant = AvogadroConstant * ElementaryCharge;

/// Speed of light (\f$c\f$) in meters per second (\f$\frac{m}{s}\f$).
const Float SpeedOfLight = 299792458;

/// Vacuum permeability (\f$\mu_{0}\f$).
const Float VacuumPermeability = 1.2566370614e-6;

/// Vacuum permittivity (\f$\epsilon_{0}\f$).
///
/// Equal to \f$\frac{1}{\mu_{0} c^{2}}\f$.
const Float VacuumPermittivity = 1.0 / (VacuumPermeability * SpeedOfLight * SpeedOfLight);

/// Fine structure constant (\f$\alpha\f$).
///
/// Equal to \f$\frac{e^{2} c \epsilon_{0}}{2 h}\f$.
const Float FineStructureConstant = (ElementaryCharge * ElementaryCharge * SpeedOfLight * VacuumPermittivity) / (2.0 * PlanckConstant);

} // end constants namespace

} // end chemkit namespace

#endif // CHEMKIT_CONSTANTS_H
