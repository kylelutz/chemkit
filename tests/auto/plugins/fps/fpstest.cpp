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

#include "fpstest.h"

#include <boost/make_shared.hpp>
#include <boost/range/algorithm.hpp>

#include <chemkit/molecule.h>
#include <chemkit/moleculefile.h>
#include <chemkit/moleculefileformat.h>

void FpsTest::initTestCase()
{
    // verify that the fps plugin registered itself correctly
    QVERIFY(boost::count(chemkit::MoleculeFileFormat::formats(), "fps") == 1);
}

void FpsTest::write()
{
    boost::shared_ptr<chemkit::Molecule> molecule =
        boost::make_shared<chemkit::Molecule>("CCO", "smiles");

    chemkit::MoleculeFile file;
    file.addMolecule(molecule);

    bool ok = file.setFormat("fps");
    QVERIFY(ok);
    QVERIFY(file.format() != 0);
    file.format()->setOption("fingerprint", "fp2");

    std::ostringstream output;
    file.write(output);

    std::string outputData = output.str();

    std::vector<std::string> lines;
    boost::split(lines,
                 outputData,
                 boost::is_any_of("\n"),
                 boost::token_compress_on);
    QCOMPARE(lines.size(), size_t(7));
    QCOMPARE(lines[0].c_str(), "#FPS1");
    QCOMPARE(lines[1].c_str(), "#num_bits=1021");
    QVERIFY(boost::starts_with(lines[2], "#type="));
    QVERIFY(boost::starts_with(lines[3], "#software=chemkit/"));
    QVERIFY(boost::starts_with(lines[4], "#date="));

    std::vector<std::string> fingerprintLine;
    boost::split(fingerprintLine,
                 lines[5],
                 boost::is_any_of("\t"),
                 boost::token_compress_on);
    QCOMPARE(fingerprintLine.size(), size_t(2));
    std::string fingerprint = fingerprintLine[0];
    QCOMPARE(fingerprint.c_str(), "0000000000000000000000000000000000000000000"
                                  "0000000000000000000000000000000000000000000"
                                  "0000000000000000000000000000000000000000000"
                                  "8000800000000000000000000000000000000400000"
                                  "0000000000000000000000000000000000000000000"
                                  "00000000000000000000000000000000000000000");

    std::string identifier = fingerprintLine[1];
    QCOMPARE(identifier.c_str(), "C2H6O");
}

QTEST_APPLESS_MAIN(FpsTest)
