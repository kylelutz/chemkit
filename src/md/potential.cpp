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

#include "potential.h"

#include <chemkit/concurrent.h>
#include <chemkit/cartesiancoordinates.h>

namespace chemkit {

// === Potential =========================================================== //
/// \class Potential potential.h chemkit/potential.h
/// \ingroup chemkit-md
/// \brief The Potential class represents a potential energy
///        expression.

// --- Construction and Destruction ---------------------------------------- //
/// Destroys the potential object.
Potential::~Potential()
{
}

// --- Properties ---------------------------------------------------------- //
/// Returns the size of the potential.
size_t Potential::size() const
{
    return 0;
}

// --- Energy -------------------------------------------------------------- //
/// Returns the potential energy of the system.
Real Potential::energy(const CartesianCoordinates *coordinates) const
{
    CHEMKIT_UNUSED(coordinates);

    return 0;
}

/// Runs the energy() method asynchronously and returns a future
/// containing the result.
///
/// \internal
boost::shared_future<Real> Potential::energyAsync(const CartesianCoordinates *coordinates) const
{
    return chemkit::concurrent::run(boost::bind(&Potential::energy, this, coordinates));
}

/// Returns the gradient of the potential energy of the system with
/// respect to \p coordinates.
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
std::vector<Vector3> Potential::gradient(const CartesianCoordinates *coordinates) const
{
   return numericalGradient(coordinates);
}

/// Returns the gradient of the potential energy of the system with
/// respect to \p coordinates. The gradient is calculated
/// numerically.
///
/// \see Potential::gradient()
std::vector<Vector3> Potential::numericalGradient(const CartesianCoordinates *coordinates) const
{
    std::vector<Vector3> gradient(coordinates->size());

    CartesianCoordinates writeableCoordinates = *coordinates;

    for(size_t i = 0; i < coordinates->size(); i++){
        const Point3 &position = coordinates->position(i);

        // initial energy
        Real eI = energy(&writeableCoordinates);
        Real epsilon = 1.0e-10;

        writeableCoordinates.setPosition(i, position + Vector3(epsilon, 0, 0));
        Real eF_x = energy(&writeableCoordinates);

        writeableCoordinates.setPosition(i, position + Vector3(0, epsilon, 0));
        Real eF_y = energy(&writeableCoordinates);

        writeableCoordinates.setPosition(i, position + Vector3(0, 0, epsilon));
        Real eF_z = energy(&writeableCoordinates);

        // restore initial position
        writeableCoordinates.setPosition(i, coordinates->position(i));

        Real dx = (eF_x - eI) / epsilon;
        Real dy = (eF_y - eI) / epsilon;
        Real dz = (eF_z - eI) / epsilon;

        gradient[i] = Vector3(dx, dy, dz);
    }

    return gradient;
}

/// Returns the root-mean-square gradient.
Real Potential::rmsg(const CartesianCoordinates *coordinates) const
{
    if(!size()){
        return 0;
    }

    Real sum = 0;

    std::vector<Vector3> gradient = this->gradient(coordinates);

    for(size_t i = 0; i < gradient.size(); i++){
        sum += gradient[i].squaredNorm();
    }

    return sqrt(sum / (3.0 * size()));
}

} // end chemkit namespace
