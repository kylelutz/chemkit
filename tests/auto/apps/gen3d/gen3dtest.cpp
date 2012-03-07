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

#include "gen3dtest.h"

#include <chemkit/atom.h>
#include <chemkit/molecule.h>
#include <chemkit/moleculefile.h>

const QString gen3dApplication = "../../../../bin/chemkit-gen3d";

void Gen3dTest::benzene()
{
    // create temporary output file
    QTemporaryFile output("XXXXXX.mol");
    output.open();

    // set arguments
    QStringList arguments;
    arguments.append("c1ccccc1");
    arguments.append(output.fileName());

    // run gen3d
    QProcess process;
    process.start(gen3dApplication, arguments);
    process.waitForFinished();

    // read output file
    QByteArray outputFileName = output.fileName().toAscii();
    chemkit::MoleculeFile file(outputFileName.constData());
    bool ok = file.read();
    if(!ok)
        qDebug() << file.errorString().c_str();
    QVERIFY(ok);

    // check formula
    boost::shared_ptr<chemkit::Molecule> molecule = file.molecule();
    QVERIFY(molecule != 0);
    QCOMPARE(molecule->formula(), std::string("C6H6"));

    // ensure the molecule has one ring
    QCOMPARE(molecule->ringCount(), size_t(1));

    // ensure center point is at (0, 0, 0)
    QVERIFY(molecule->center().isZero(0.01));

    // ensure each atom position is not at (0, 0, 0)
    foreach(const chemkit::Atom *atom, molecule->atoms()){
        QVERIFY(!atom->position().isZero(0.01));
    }
}

QTEST_APPLESS_MAIN(Gen3dTest)
