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

#ifndef CHEMKIT_AMINOACID_H
#define CHEMKIT_AMINOACID_H

#include "chemkit.h"

#include "residue.h"
#include "vector3.h"

namespace chemkit {

class AminoAcidPrivate;

class CHEMKIT_EXPORT AminoAcid : public Residue
{
    public:
        enum Conformation {
            Coil,
            AlphaHelix,
            BetaSheet
        };

        enum AminoAcidType {
            Alanine,
            Arganine,
            Asparagine,
            AsparticAcid,
            Aspartate = AsparticAcid,
            Cystenine,
            Glutamine,
            GlutamicAcid,
            Glutamate = GlutamicAcid,
            Glycine,
            Histadine,
            Isoleucine,
            Leucine,
            Lysine,
            Methionine,
            Phenylalanine,
            Proline,
            Serine,
            Threonine,
            Tryptophan,
            Tyrosine,
            Valine,
            UnspecifiedType
        };

        // construction and destruction
        AminoAcid(Molecule *molecule);
        ~AminoAcid();

        // properties
        void setType(AminoAcidType type);
        void setType(const QString &letterOrSymbol);
        AminoAcidType type() const;
        QString name() const;
        QString symbol() const;
        QString letter() const;
        void setConformation(Conformation conformation);
        Conformation conformation() const;

        // structure
        void setAlphaCarbon(Atom *atom);
        Atom* alphaCarbon() const;
        void setCarbonylCarbon(Atom *atom);
        Atom* carbonylCarbon() const;
        void setCarbonylOxygen(Atom *atom);
        Atom* carbonylOxygen() const;
        void setAminoNitrogen(Atom *atom);
        Atom* aminoNitrogen() const;

        // geometry
        Vector3 peptidePlaneNormal() const;

    private:
        AminoAcidPrivate* const d;
};

} // end chemkit namespace

#endif // CHEMKIT_AMINOACID_H
