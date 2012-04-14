/******************************************************************************
**
** Copyright (C) 2009-2012 Kyle Lutz <kyle.r.lutz@gmail.com>
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

#include "integrator.h"

#include <chemkit/cartesiancoordinates.h>

#include "potential.h"

namespace chemkit {

// === IntegratorPrivate =================================================== //
class IntegratorPrivate
{
public:
    boost::shared_ptr<Potential> potential;
    CartesianCoordinates coordinates;
};

// === Integrator ========================================================== //
/// \class Integrator integrator.h chemkit/integrator.h
/// \ingroup chemkit-md
/// \brief The Integrator class represents an integrator.

// --- Construction and Destruction ---------------------------------------- //
/// Creates a new integrator.
Integrator::Integrator()
    : d(new IntegratorPrivate)
{
}

/// Destroys the integrator object.
Integrator::~Integrator()
{
    delete d;
}

// --- Potential ----------------------------------------------------------- //
/// Sets the potential for the integrator to \p potential.
void Integrator::setPotential(const boost::shared_ptr<Potential> &potential)
{
    d->potential = potential;
}

/// Returns the potential for the integrator.
boost::shared_ptr<Potential> Integrator::potential() const
{
    return d->potential;
}

// --- Coordinates --------------------------------------------------------- //
/// Sets the initial coordinates to \p coordinates.
void Integrator::setCoordinates(const CartesianCoordinates *coordinates)
{
    d->coordinates = *coordinates;
}

/// Returns the current coordinates.
CartesianCoordinates* Integrator::coordinates() const
{
    return &d->coordinates;
}

// --- Energy -------------------------------------------------------------- //
/// Returns the energy of the system.
Real Integrator::energy() const
{
    if(!d->potential){
        return 0;
    }

    return d->potential->energy(&d->coordinates);
}

/// Returns the gradient of the energy.
std::vector<Vector3> Integrator::gradient() const
{
    if(!d->potential){
        return std::vector<Vector3>();
    }

    return d->potential->gradient(&d->coordinates);
}

/// Returns the root-mean-square gradient.
Real Integrator::rmsg() const
{
    if(!d->potential){
        return 0;
    }

    return d->potential->rmsg(&d->coordinates);
}

// --- Integration --------------------------------------------------------- //
/// Performs a single integration step.
void Integrator::integrate()
{
}

} // end chemkit namespace
