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

#ifndef CHEMKIT_POLYMERCHAIN_H
#define CHEMKIT_POLYMERCHAIN_H

#include "chemkit.h"

#include <string>

#include <QtCore>

namespace chemkit {

class Polymer;
class Residue;
class PolymerChainPrivate;

class CHEMKIT_EXPORT PolymerChain
{
    public:
        // properties
        void setName(const std::string &name);
        std::string name() const;
        int size() const;
        bool isEmpty() const;
        Polymer* polymer() const;

        // structure
        void addResidue(Residue *residue);
        void appendResidue(Residue *residue);
        void prependResidue(Residue *residue);
        void insertResidue(int index, Residue *residue);
        bool removeResidue(Residue *residue);
        bool deleteResidue(Residue *residue);
        Residue* residue(int index) const;
        QList<Residue *> residues() const;
        int residueCount() const;
        int indexOf(const Residue *residue) const;
        std::string sequenceString() const;
        int sequenceNumber(const Residue *residue) const;

    private:
        // construction and destruction
        PolymerChain(Polymer *polymer);
        ~PolymerChain();

        Q_DISABLE_COPY(PolymerChain)

        friend class Polymer;

    private:
        PolymerChainPrivate* const d;
};

} // end chemkit namespace

#endif // CHEMKIT_POLYMERCHAIN_H
