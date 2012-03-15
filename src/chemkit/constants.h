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
const Real Pi = 3.14159265358979323846264338327;

/// Tau (\f$\tau\f$). Twice the number pi (\f$2\pi\f$). The ratio
/// of the circumference of a circle to its radius.
const Real Tau = 2.0 * Pi;

/// Conversion factor from degrees to radians.
///
/// Equal to \f$\frac{\pi}{180}\f$.
const Real DegreesToRadians = Pi / 180.0;

/// Conversion factor from radians to degrees.
///
/// Equal to \f$\frac{180}{\pi}\f$.
const Real RadiansToDegrees = 180.0 / Pi;

/// Conversion factor from calories to joules.
const Real CaloriesToJoules = 4.184;

/// Conversion factor from joules to calories.
const Real JoulesToCalories = 1.0 / CaloriesToJoules;

/// Avogadro constant (\f$N_{A}\f$).
const Real AvogadroConstant = 6.02214179e23;

/// Boltzmann constant (\f$k_{B}\f$)
const Real BoltzmannConstant = 1.3806503e-23;

/// Gas constant.
///
/// Equal to \f$N_{A} k_{B}\f$.
const Real GasConstant = AvogadroConstant * BoltzmannConstant;

/// Planck constant (\f$h\f$).
const Real PlanckConstant = 6.62606896e-34;

/// Reduced Planck constant (\f$\hbar\f$).
const Real ReducedPlanckConstant = PlanckConstant / (2.0 * Pi);

/// Mass of a proton (\f$m_{p}\f$).
const Real ProtonMass = 1.672621637e-27;

/// Mass of a electron (\f$m_{e}\f$).
const Real ElectronMass = 9.10938215e-31;

/// Elementary charge (\f$e\f$).
const Real ElementaryCharge = 1.602176487e-19;

/// Charge of a proton (\f$+e\f$).
const Real ProtonCharge = ElementaryCharge;

/// Charge of a electron (\f$-e\f$).
const Real ElectronCharge = -ElementaryCharge;

/// Faraday constant (\f$F\f$).
///
/// Equal to \f$e N_{A}\f$.
const Real FaradayConstant = AvogadroConstant * ElementaryCharge;

/// Speed of light (\f$c\f$) in meters per second (\f$\frac{m}{s}\f$).
const Real SpeedOfLight = 299792458;

/// Vacuum permeability (\f$\mu_{0}\f$).
const Real VacuumPermeability = 1.2566370614e-6;

/// Vacuum permittivity (\f$\epsilon_{0}\f$).
///
/// Equal to \f$\frac{1}{\mu_{0} c^{2}}\f$.
const Real VacuumPermittivity = 1.0 / (VacuumPermeability * SpeedOfLight * SpeedOfLight);

/// Fine structure constant (\f$\alpha\f$).
///
/// Equal to \f$\frac{e^{2} c \epsilon_{0}}{2 h}\f$.
const Real FineStructureConstant = (ElementaryCharge * ElementaryCharge * SpeedOfLight * VacuumPermittivity) / (2.0 * PlanckConstant);

/// Conversion factor between Bohr units and Angstroms.
const Real BohrToAnstroms = 0.52918;

} // end constants namespace

} // end chemkit namespace

#endif // CHEMKIT_CONSTANTS_H
