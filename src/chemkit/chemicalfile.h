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

#ifndef CHEMKIT_CHEMICALFILE_H
#define CHEMKIT_CHEMICALFILE_H

#include "chemkit.h"

#include <string>
#include <vector>

#include <QtCore>

namespace chemkit {

class Molecule;
class ChemicalFileFormat;
class ChemicalFilePrivate;

class CHEMKIT_EXPORT ChemicalFile
{
    public:
        // construction and destruction
        ChemicalFile();
        ChemicalFile(const std::string &fileName);
        ~ChemicalFile();

        // properties
        void setFileName(const std::string &fileName);
        std::string fileName() const;
        void setFormat(ChemicalFileFormat *format);
        bool setFormat(const std::string &name);
        ChemicalFileFormat* format() const;
        std::string formatName() const;
        int size() const;
        bool isEmpty() const;

        // file contents
        void addMolecule(Molecule *molecule);
        bool removeMolecule(Molecule *molecule);
        bool deleteMolecule(Molecule *molecule);
        QList<Molecule *> molecules() const;
        int moleculeCount() const;
        Molecule* molecule(int index = 0) const;
        bool contains(const Molecule *molecule) const;
        void clear();

        // file data
        void setFileData(const QString &name, const QVariant &value);
        QVariant fileData(const QString &name) const;
        QHash<QString, QVariant> fileData() const;
        void setMoleculeData(const Molecule *molecule, const QString &name, const QVariant &value);
        QVariant moleculeData(const Molecule *molecule, const QString &name) const;
        QHash<QString, QVariant> moleculeData(const Molecule *molecule) const;

        // input and output
        bool read();
        bool read(const std::string &fileName);
        bool read(const std::string &fileName, const std::string &format);
        bool read(QIODevice *iodev, const std::string &format);
        bool write();
        bool write(const std::string &fileName);
        bool write(const std::string &fileName, const std::string &format);
        bool write(QIODevice *iodev);
        bool write(QIODevice *iodev, const std::string &format);

        // error handling
        std::string errorString() const;

        // static methods
        static std::vector<std::string> formats();
        static Molecule* quickRead(const std::string &fileName);
        static void quickWrite(const Molecule *molecule, const std::string &fileName);

    private:
        void setErrorString(const std::string &error);

    private:
        ChemicalFilePrivate* const d;
};

} // end chemkit namespace

#endif // CHEMKIT_CHEMICALFILE_H
