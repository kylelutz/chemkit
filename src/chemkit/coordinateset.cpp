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

#include "coordinateset.h"

#include <boost/scoped_ptr.hpp>

#include "diagramcoordinates.h"
#include "internalcoordinates.h"
#include "cartesiancoordinates.h"

namespace chemkit {

// === CoordinateSet ======================================================= //
/// \class CoordinateSet coordinateset.h chemkit/coordinateset.h
/// \ingroup chemkit
/// \brief The CoordinateSet class contains a set of coordinates.

// --- Construction and Destruction ---------------------------------------- //
/// Create a new coordinate set with no coordinates.
CoordinateSet::CoordinateSet()
{
    m_type = None;
    m_cartesianCordinates = 0;
}

/// Creates a new coordinate set with \p coordinates.
CoordinateSet::CoordinateSet(CartesianCoordinates *coordinates)
{
    m_type = Cartesian;
    m_cartesianCordinates = coordinates;
}

/// Creates a new coordinate set with \p coordinates.
CoordinateSet::CoordinateSet(InternalCoordinates *coordinates)
{
    m_type = Internal;
    m_internalCoordinates = coordinates;
}

/// Creates a new coordinate set with \p coordinates.
CoordinateSet::CoordinateSet(DiagramCoordinates *coordinates)
{
    m_type = Diagram;
    m_diagramCoordinates = coordinates;
}

/// Creates a new coordinate set as a copy of \p other.
CoordinateSet::CoordinateSet(const CoordinateSet &other)
{
    m_type = None;

    if(other.type() == Cartesian){
        setCoordinates(new CartesianCoordinates(*other.cartesianCoordinates()));
    }
    else if(other.type() == Internal){
        setCoordinates(new InternalCoordinates(*other.internalCoordinates()));
    }
    else if(other.type() == Diagram){
        setCoordinates(new DiagramCoordinates(*other.diagramCoordinates()));
    }
}

/// Destroys the coordinate set object.
CoordinateSet::~CoordinateSet()
{
    clear();
}

// --- Properties ---------------------------------------------------------- //
/// Returns the type of coordinates that the coordinate set contains.
CoordinateSet::Type CoordinateSet::type() const
{
    return m_type;
}

/// Returns the size of the coordinates stored in the coordinate
/// set.
size_t CoordinateSet::size() const
{
    switch(m_type){
        case Cartesian: return m_cartesianCordinates->size();
        case Internal: return m_internalCoordinates->size();
        case Diagram: return m_diagramCoordinates->size();
        default: return 0;
    }
}

/// Returns \c true if the coordinates stored are empty.
bool CoordinateSet::isEmpty() const
{
    return size() == 0;
}

/// Sets the type to \c Cartesian and the coordinates to
/// \p coordinates.
void CoordinateSet::setCoordinates(CartesianCoordinates *coordinates)
{
    clear();

    m_type = Cartesian;
    m_cartesianCordinates = coordinates;
}

/// Sets the type to \c Internal and the coordinates to
/// \p coordinates.
void CoordinateSet::setCoordinates(InternalCoordinates *coordinates)
{
    clear();

    m_type = Internal;
    m_internalCoordinates = coordinates;
}

/// Sets the type to \c Diagram and the coordinates to
/// \p coordinates.
void CoordinateSet::setCoordinates(DiagramCoordinates *coordinates)
{
    clear();

    m_type = Diagram;
    m_diagramCoordinates = coordinates;
}

/// Returns the cartesian coorinates stored in the coordinate
/// set. Returns \c 0 if the coordinate set does not contain
/// cartesian coordinates.
CartesianCoordinates* CoordinateSet::cartesianCoordinates() const
{
    if(m_type == Cartesian){
        return m_cartesianCordinates;
    }

    return 0;
}

/// Returns the internal coorinates stored in the coordinate
/// set. Returns \c 0 if the coordinate set does not contain
/// internal coordinates.
InternalCoordinates* CoordinateSet::internalCoordinates() const
{
    if(m_type == Internal){
        return m_internalCoordinates;
    }

    return 0;
}

/// Returns the diagram coorinates stored in the coordinate
/// set. Returns \c 0 if the coordinate set does not contain
/// diagram coordinates.
DiagramCoordinates* CoordinateSet::diagramCoordinates() const
{
    if(m_type == Diagram){
        return m_diagramCoordinates;
    }

    return 0;
}

/// Clears the coordinates stored in the coordinate set.
void CoordinateSet::clear()
{
    if(m_type == Cartesian){
        delete m_cartesianCordinates;
        m_cartesianCordinates = 0;
    }
    else if(m_type == Internal){
        delete m_internalCoordinates;
        m_internalCoordinates = 0;
    }
    else if(m_type == Diagram){
        delete m_diagramCoordinates;
        m_diagramCoordinates = 0;
    }

    m_type = None;
}

// --- Position ------------------------------------------------------------ //
/// Returns the 3D cartesian position of the point at \p index.
Point3 CoordinateSet::position(size_t index) const
{
    if(m_type == Cartesian){
        return m_cartesianCordinates->position(index);
    }
    else if(m_type == Internal){
        boost::scoped_ptr<CartesianCoordinates>
            coordinates(m_internalCoordinates->toCartesianCoordinates());

        return coordinates->position(index);
    }
    else if(m_type == Diagram){
        Point2f point2 = m_diagramCoordinates->position(index);

        return Point3(point2.x(), point2.y(), 0);
    }

    return Point3();
}

// --- Operators ----------------------------------------------------------- //
CoordinateSet& CoordinateSet::operator=(const CoordinateSet &other)
{
    if(this != &other){
        if(other.type() == None){
            clear();
        }
        else if(other.type() == Cartesian){
            setCoordinates(new CartesianCoordinates(*other.cartesianCoordinates()));
        }
        else if(other.type() == Internal){
            setCoordinates(new InternalCoordinates(*other.internalCoordinates()));
        }
        else if(other.type() == Diagram){
            setCoordinates(new DiagramCoordinates(*other.diagramCoordinates()));
        }
    }

    return *this;
}

} // end chemkit namespace
