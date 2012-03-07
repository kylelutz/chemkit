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

#include "gromacstest.h"

#include <boost/range/algorithm.hpp>

#include <chemkit/topology.h>
#include <chemkit/topologyfile.h>
#include <chemkit/topologyfileformat.h>

const std::string dataPath = "../../../data/";

void GromacsTest::initTestCase()
{
    // verify that the gromacs plugin registered itself correctly
    QVERIFY(boost::count(chemkit::TopologyFileFormat::formats(), "gro") == 1);
    QVERIFY(boost::count(chemkit::TopologyFileFormat::formats(), "top") == 1);
}

void GromacsTest::spc216()
{
    chemkit::TopologyFile file(dataPath + "spc216.gro");
    bool ok = file.read();
    if(!ok)
        qDebug() << file.errorString().c_str();
    QVERIFY(ok);

    boost::shared_ptr<chemkit::Topology> topology = file.topology();
    QVERIFY(topology != 0);
    QCOMPARE(topology->size(), size_t(648));

    for(size_t i = 0; i < 648; i += 3){
        QCOMPARE(topology->type(i+0), std::string("OW"));
        QCOMPARE(topology->type(i+1), std::string("HW1"));
        QCOMPARE(topology->type(i+2), std::string("HW2"));
    }
}

void GromacsTest::ubiquitin()
{
    chemkit::TopologyFile file(dataPath + "1UBQ.top");
    bool ok = file.read();
    if(!ok)
        qDebug() << file.errorString().c_str();
    QVERIFY(ok);

    boost::shared_ptr<chemkit::Topology> topology = file.topology();
    QVERIFY(topology != 0);
    QCOMPARE(topology->size(), size_t(1231));
}

QTEST_APPLESS_MAIN(GromacsTest)
