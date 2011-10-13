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

#ifndef CHEMKIT_INTERNALCOORDINATES_H
#define CHEMKIT_INTERNALCOORDINATES_H

#include "chemkit.h"

#include <vector>

#include "point3.h"

namespace chemkit {

class Coordinates;
class InternalCoordinatesPrivate;

class CHEMKIT_EXPORT InternalCoordinates
{
    public:
        // construction and destruction
        InternalCoordinates();
        InternalCoordinates(int size);
        InternalCoordinates(const InternalCoordinates &coordinates);
        ~InternalCoordinates();

        // properties
        int size() const;

        // coordinates
        void setCoordinates(int row, Real r, Real theta = 0, Real phi = 0);
        void setCoordinatesRadians(int row, Real r, Real theta = 0, Real phi = 0);
        std::vector<Real> coordinates(int row) const;
        std::vector<Real> coordinatesRadians(int row) const;
        void setConnections(int row, int a, int b = 0, int c = 0);
        std::vector<int> connections(int row) const;

        // conversions
        Coordinates* toCartesianCoordinates() const;

        // operators
        InternalCoordinates& operator=(const InternalCoordinates &coordinates);

    private:
        InternalCoordinatesPrivate* const d;
};

} // end chemkit namespace

#endif // CHEMKIT_INTERNALCOORDINATES_H
