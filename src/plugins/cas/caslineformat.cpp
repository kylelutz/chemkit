/******************************************************************************
**
** Copyright (C) 2009-2012 Kyle Lutz <kyle.r.lutz@gmail.com>
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

#include "caslineformat.h"

#include <QtCore>

#include <chemkit/pubchem.h>
#include <chemkit/molecule.h>

CasLineFormat::CasLineFormat()
    : chemkit::LineFormat("cas")
{
}

chemkit::Molecule* CasLineFormat::read(const std::string &formula)
{
    // ensure a QApplication object has been created
    int argc = 0;
    QCoreApplication app(argc, 0);
    Q_UNUSED(app);

    // lookup the compound id for the cas formula
    chemkit::PubChem pubchem;
    QStringList ids = pubchem.search(formula.c_str());

    if(ids.isEmpty()){
        setErrorString("Failed to lookup CAS number");
        return 0;
    }

    // download molecule from its compound id
    boost::shared_ptr<chemkit::Molecule> downloadedMolecule = pubchem.downloadMolecule(ids.first());

    if(!downloadedMolecule){
        setErrorString("Failed to download molecule");
        return 0;
    }

    return new chemkit::Molecule(*downloadedMolecule);
}
