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

#ifndef CHEMKIT_SCALARFIELD_H
#define CHEMKIT_SCALARFIELD_H

#include "chemkit.h"

#include <vector>

#include "point3.h"
#include "vector3.h"

namespace chemkit {

class ScalarFieldPrivate;

class CHEMKIT_EXPORT ScalarField
{
    public:
        // construction and destruction
        ScalarField();
        ScalarField(const std::vector<int> &dimensions, const std::vector<Real> &cellLengths, const std::vector<Real> &data);
        ~ScalarField();

        // properties
        int width() const;
        int height() const;
        int depth() const;
        int size() const;
        std::vector<int> dimensions() const;
        Real cellWidth() const;
        Real cellHeight() const;
        Real cellDepth() const;
        std::vector<Real> cellDimensions() const;
        void setOrigin(const Point3 &origin);
        Point3 origin() const;
        std::vector<Real> data() const;

        // values
        void setValue(int i, int j, int k, Real value);
        Real value(int i, int j, int k) const;
        Real value(const Point3 &position) const;
        Point3 position(int i, int j, int k) const;
        Vector3 gradient(int i, int j, int k) const;
        Vector3 gradient(const Point3 &position) const;

    private:
        ScalarFieldPrivate* const d;
};

} // end chemkit namespace

#endif // CHEMKIT_SCALARFIELD_H
