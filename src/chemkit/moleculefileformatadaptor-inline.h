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

#ifndef CHEMKIT_MOLECULEFILEFORMATADAPTOR_INLINE_H
#define CHEMKIT_MOLECULEFILEFORMATADAPTOR_INLINE_H

#include "moleculefileformatadaptor.h"

#include <boost/foreach.hpp>

#include "polymer.h"
#include "lineformat.h"
#include "polymerfile.h"
#include "polymerfileformat.h"

namespace chemkit {

// === MoleculeFileFormatAdaptor<LineFormat> =============================== //
inline MoleculeFileFormatAdaptor<LineFormat>::MoleculeFileFormatAdaptor(LineFormat *format)
    : MoleculeFileFormat(format->name())
{
    m_format = format;
}

inline MoleculeFileFormatAdaptor<LineFormat>::~MoleculeFileFormatAdaptor()
{
    delete m_format;
}

inline bool MoleculeFileFormatAdaptor<LineFormat>::read(QIODevice *iodev, MoleculeFile *file)
{
    iodev->setTextModeEnabled(true);

    while(!iodev->atEnd()){
        QString line = iodev->readLine().simplified();

        QStringList lineItems = line.split(' ', QString::SkipEmptyParts);
        if(lineItems.isEmpty()){
            continue;
        }

        std::string formula = lineItems[0].toStdString();
        Molecule *molecule = m_format->read(formula);
        if(!molecule){
            continue;
        }

        if(lineItems.size() >= 2){
            std::string name = line.mid(formula.length()).trimmed().toStdString();
            molecule->setName(name);
        }

        file->addMolecule(molecule);
    }

    return true;
}

inline bool MoleculeFileFormatAdaptor<LineFormat>::write(const MoleculeFile *file, QIODevice *iodev)
{
    iodev->setTextModeEnabled(true);

    BOOST_FOREACH(const Molecule *molecule, file->molecules()){
        std::string formula = m_format->write(molecule);
        bool ok = iodev->write(formula.c_str());
        if(!ok){
            continue;
        }

        if(!molecule->name().empty()){
            iodev->write(" ");
            iodev->write(molecule->name().c_str());
        }

        iodev->write("\n");
    }

    return true;
}

// === MoleculeFileFormatAdaptor<PolymerFileFormat> ======================= //
inline MoleculeFileFormatAdaptor<PolymerFileFormat>::MoleculeFileFormatAdaptor(PolymerFileFormat *format)
    : MoleculeFileFormat(format->name())
{
    m_format = format;
}

inline MoleculeFileFormatAdaptor<PolymerFileFormat>::~MoleculeFileFormatAdaptor()
{
    delete m_format;
}

inline bool MoleculeFileFormatAdaptor<PolymerFileFormat>::read(QIODevice *iodev, MoleculeFile *file)
{
    PolymerFile polymerFile;
    bool ok = polymerFile.read(iodev, m_format->name());
    if(!ok){
        setErrorString(polymerFile.errorString());
        return false;
    }

    Q_FOREACH(Polymer *polymer, polymerFile.polymers()){
        // remove polymer from the polymer file
        polymerFile.removePolymer(polymer);

        // add polymer to the molecule file
        file->addMolecule(polymer);
    }

    return true;
}

} // end chemkit namespace

#endif // CHEMKIT_MOLECULEFILEFORMATADAPTOR_INLINE_H
