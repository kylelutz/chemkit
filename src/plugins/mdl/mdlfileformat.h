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

#ifndef MDLFILEFORMAT_H
#define MDLFILEFORMAT_H

#include <QtCore>

#include <chemkit/molecule.h>
#include <chemkit/moleculefileformat.h>

class MdlFileFormat : public chemkit::MoleculeFileFormat
{
public:
    // construction and destruction
    MdlFileFormat(const std::string &name);
    ~MdlFileFormat();

    // input and output
    bool read(std::istream &input, chemkit::MoleculeFile *file);
    bool read(QIODevice *iodev, chemkit::MoleculeFile *file);
    bool write(const chemkit::MoleculeFile *file, std::ostream &output);
    bool write(const chemkit::MoleculeFile *file, QIODevice *iodev);

private:
    bool readMolFile(QIODevice *iodev, chemkit::MoleculeFile *file);
    bool readSdfFile(QIODevice *iodev, chemkit::MoleculeFile *file);
    bool readAtomBlock(QIODevice *iodev, chemkit::Molecule *molecule, int atomCount);
    bool readBondBlock(QIODevice *iodev, chemkit::Molecule *molecule, int bondCount);
    bool readPropertyBlock(QIODevice *iodev, chemkit::Molecule *molecule);
    bool readDataBlock(QIODevice *iodev, chemkit::Molecule *molecule, chemkit::MoleculeFile *file);
    void writeMolFile(const chemkit::Molecule *molecule, QIODevice *iodev);
    void writeSdfFile(const chemkit::MoleculeFile *file, QIODevice *iodev);
    void writeAtomBlock(const chemkit::Molecule *molecule, QIODevice *iodev);
    void writeBondBlock(const chemkit::Molecule *molecule, QIODevice *iodev);
};

#endif // MDLFILEFORMAT_H
