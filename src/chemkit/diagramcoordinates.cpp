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

#include "diagramcoordinates.h"

#include "cartesiancoordinates.h"

namespace chemkit {

// === DiagramCoordinates ================================================== //
/// \class DiagramCoordinates diagramcoordinates.h chemkit/diagramcoordinates.h
/// \ingroup chemkit
/// \brief The DiagramCoordinates class contains 2D positions for
///        atoms in a molecular diagram.

// --- Construction and Destruction ---------------------------------------- //
/// Creates a new diagram coordinates object with space for \p size
/// points.
DiagramCoordinates::DiagramCoordinates(size_t size)
    : m_coordinates(size)
{
}

/// Destroys the diagram coordinates object.
DiagramCoordinates::~DiagramCoordinates()
{
}

// --- Properties ---------------------------------------------------------- //
/// Sets the number of points in the coordinates to \p size.
void DiagramCoordinates::resize(size_t size)
{
    m_coordinates.resize(size);
}

/// Returns the size of the diagram coordinates.
size_t DiagramCoordinates::size() const
{
    return m_coordinates.size();
}

/// Returns \c true if the diagram coordinates contains no points.
bool DiagramCoordinates::isEmpty() const
{
    return m_coordinates.empty();
}

// --- Positions ----------------------------------------------------------- //
/// Set the position at \p index to \p position.
void DiagramCoordinates::setPosition(size_t index, const Point2f &position)
{
    assert(index < m_coordinates.size());

    m_coordinates[index] = position;
}

/// Sets the position at \p index to (\p x, \p y).
void DiagramCoordinates::setPosition(size_t index, float x, float y)
{
    setPosition(index, Point2f(x, y));
}

/// Returns the position at \p index.
Point2f DiagramCoordinates::position(size_t index) const
{
    assert(index < m_coordinates.size());

    return m_coordinates[index];
}

/// Adds \p position to the diagram coordinates.
void DiagramCoordinates::append(const Point2f &position)
{
    m_coordinates.push_back(position);
}

/// Inserts \p position at \p index.
void DiagramCoordinates::insert(size_t index, const Point2f &position)
{
    m_coordinates.insert(m_coordinates.begin() + index, position);
}

/// Removes the point at \p index.
void DiagramCoordinates::remove(size_t index)
{
    m_coordinates.erase(m_coordinates.begin() + index);
}

// --- Conversions --------------------------------------------------------- //
/// Converts the diagram coordinates into 3D cartesian coordinates.
///
/// The ownership of the returned coordinates object is passed to
/// the caller.
CartesianCoordinates* DiagramCoordinates::toCartesianCoordinates() const
{
    CartesianCoordinates *coordinates = new CartesianCoordinates(m_coordinates.size());

    for(size_t i = 0; i < m_coordinates.size(); i++){
        const Point2f &point = m_coordinates[i];

        coordinates->setPosition(i, Point3(point.x(), point.y(), 0));
    }

    return coordinates;
}

} // end chemkit namespace
