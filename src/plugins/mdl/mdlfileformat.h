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

#ifndef MDLFILEFORMAT_H
#define MDLFILEFORMAT_H

#include <QtCore>

#include <chemkit/molecule.h>
#include <chemkit/chemicalfileformat.h>

class MdlFileFormat : public chemkit::ChemicalFileFormat
{
    public:
        // construction and destruction
        MdlFileFormat(const QString &name);
        ~MdlFileFormat();

        // input and output
        bool read(QIODevice *iodev, chemkit::ChemicalFile *file);
        bool write(const chemkit::ChemicalFile *file, QIODevice *iodev);

    private:
        bool readMolFile(QIODevice *iodev, chemkit::ChemicalFile *file);
        bool readSdfFile(QIODevice *iodev, chemkit::ChemicalFile *file);
        bool readAtomBlock(QIODevice *iodev, chemkit::Molecule *molecule, int atomCount);
        bool readBondBlock(QIODevice *iodev, chemkit::Molecule *molecule, int bondCount);
        bool readPropertyBlock(QIODevice *iodev, chemkit::Molecule *molecule);
        bool readDataBlock(QIODevice *iodev, const chemkit::Molecule *molecule, chemkit::ChemicalFile *file);
        void writeMolFile(const chemkit::Molecule *molecule, QIODevice *iodev);
        void writeSdfFile(const chemkit::ChemicalFile *file, QIODevice *iodev);
        void writeAtomBlock(const chemkit::Molecule *molecule, QIODevice *iodev);
        void writeBondBlock(const chemkit::Molecule *molecule, QIODevice *iodev);
};

#endif // MDLFILEFORMAT_H
