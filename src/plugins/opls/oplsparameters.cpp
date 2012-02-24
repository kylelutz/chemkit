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

#include "oplsparameters.h"

#include <fstream>

#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>

#include <chemkit/foreach.h>

// --- Construction and Destruction ---------------------------------------- //
OplsParameters::OplsParameters(const std::string &fileName)
    : m_fileName(fileName)
{
    read(fileName);
}

OplsParameters::~OplsParameters()
{
}

// --- Properties ---------------------------------------------------------- //
void OplsParameters::setFileName(const std::string &fileName)
{
    m_fileName = fileName;
}

std::string OplsParameters::fileName() const
{
    return m_fileName;
}

// --- Parameters ---------------------------------------------------------- //
int OplsParameters::atomClass(int type) const
{
    return m_typeToClass[type];
}

std::string OplsParameters::atomName(int type) const
{
    return m_typeToName[type];
}

chemkit::Real OplsParameters::partialCharge(int type) const
{
    return m_typeToCharge[type];
}

const OplsBondStrechParameters* OplsParameters::bondStrechParameters(int a, int b) const
{
    a = atomClass(a);
    b = atomClass(b);

    if(a > b){
        std::swap(a, b);
    }

    foreach(const OplsBondStrechParameters &p, m_bondStrechParameters){
        if(p.typeA == a && p.typeB == b){
            return &p;
        }
    }

    return 0;
}

const OplsAngleBendParameters* OplsParameters::angleBendParameters(int a, int b, int c) const
{
    a = atomClass(a);
    b = atomClass(b);
    c = atomClass(c);

    if(a > c){
        std::swap(a, c);
    }

    foreach(const OplsAngleBendParameters &p, m_angleBendParameters){
        if(p.typeA == a && p.typeB == b && p.typeC == c){
            return &p;
        }
    }

    return 0;
}

const OplsTorsionParameters* OplsParameters::torsionParameters(int a, int b, int c, int d) const
{
    a = atomClass(a);
    b = atomClass(b);
    c = atomClass(c);
    d = atomClass(d);

    if(b > c){
        std::swap(b, c);
    }

    if(a > d){
        std::swap(a, d);
    }

    foreach(const OplsTorsionParameters &p, m_torsionParameters){
        if(p.typeA == a && p.typeB == b && p.typeC == c && p.typeD == d){
            return &p;
        }
    }

    return 0;
}

const OplsVanDerWaalsParameters* OplsParameters::vanDerWaalsParameters(int type) const
{
    return &m_vanDerWaalsParameters[type];
}

// --- Internal Methods ---------------------------------------------------- //
bool OplsParameters::read(const std::string &fileName)
{
    std::ifstream file(fileName.c_str());
    if(!file.is_open()){
        return false;
    }

    while(!file.eof()){
        std::string line;
        std::getline(file, line);

        // atom parameters
        if(boost::starts_with(line, "atom")){
            std::vector<std::string> lineItems;
            boost::split(lineItems,
                         line,
                         boost::is_any_of(" \t"),
                         boost::token_compress_on);

            if(lineItems.size() < 4){
                continue;
            }

            int type = boost::lexical_cast<int>(lineItems[1]);
            int klass = boost::lexical_cast<int>(lineItems[2]);
            std::string name = lineItems[3];

            if(m_typeToClass.size() < static_cast<size_t>(type + 1))
                m_typeToClass.resize(type + 1);
            m_typeToClass[type] = klass;

            if(m_typeToName.size() < static_cast<size_t>(type + 1))
                m_typeToName.resize(type + 1);
            m_typeToName[type] = name;
        }
        // bond parameters
        else if(boost::starts_with(line, "bond")){
            std::vector<std::string> lineItems;
            boost::split(lineItems,
                         line,
                         boost::is_any_of(" \t"),
                         boost::token_compress_on);

            if(lineItems.size() < 5){
                continue;
            }

            OplsBondStrechParameters p;
            p.typeA = boost::lexical_cast<int>(lineItems[1]);
            p.typeB = boost::lexical_cast<int>(lineItems[2]);
            p.kb = boost::lexical_cast<chemkit::Real>(lineItems[3]);
            p.r0 = boost::lexical_cast<chemkit::Real>(lineItems[4]);
            m_bondStrechParameters.push_back(p);
        }
        // angle parameters
        else if(boost::starts_with(line, "angle")){
            std::vector<std::string> lineItems;
            boost::split(lineItems,
                         line,
                         boost::is_any_of(" \t"),
                         boost::token_compress_on);

            if(lineItems.size() < 6){
                continue;
            }

            OplsAngleBendParameters p;
            p.typeA = boost::lexical_cast<int>(lineItems[1]);
            p.typeB = boost::lexical_cast<int>(lineItems[2]);
            p.typeC = boost::lexical_cast<int>(lineItems[3]);
            p.ka = boost::lexical_cast<chemkit::Real>(lineItems[4]);
            p.theta0 = boost::lexical_cast<chemkit::Real>(lineItems[5]);
            m_angleBendParameters.push_back(p);
        }
        // torsion parameters
        else if(boost::starts_with(line, "torsion")){
            std::vector<std::string> lineItems;
            boost::split(lineItems,
                         line,
                         boost::is_any_of(" \t"),
                         boost::token_compress_on);

            if(lineItems.size() < 14){
                continue;
            }

            OplsTorsionParameters p;
            p.typeA = boost::lexical_cast<int>(lineItems[1]);
            p.typeB = boost::lexical_cast<int>(lineItems[2]);
            p.typeC = boost::lexical_cast<int>(lineItems[3]);
            p.typeD = boost::lexical_cast<int>(lineItems[4]);
            p.v1 = boost::lexical_cast<chemkit::Real>(lineItems[5]);
            p.v2 = boost::lexical_cast<chemkit::Real>(lineItems[8]);
            p.v3 = boost::lexical_cast<chemkit::Real>(lineItems[11]);
            m_torsionParameters.push_back(p);
        }
        // van der waals parameters
        else if(boost::starts_with(line, "vdw")){
            std::vector<std::string> lineItems;
            boost::split(lineItems,
                         line,
                         boost::is_any_of(" \t"),
                         boost::token_compress_on);

            if(lineItems.size() < 4){
                continue;
            }

            int type = boost::lexical_cast<int>(lineItems[1]);
            if(m_vanDerWaalsParameters.size() < static_cast<size_t>(type + 1))
                m_vanDerWaalsParameters.resize(type + 1);

            OplsVanDerWaalsParameters p;
            p.sigma = boost::lexical_cast<chemkit::Real>(lineItems[2]);
            p.epsilon = boost::lexical_cast<chemkit::Real>(lineItems[3]);
            m_vanDerWaalsParameters[type] = p;
        }
        else if(boost::starts_with(line, "charge")){
            std::vector<std::string> lineItems;
            boost::split(lineItems,
                         line,
                         boost::is_any_of(" \t"),
                         boost::token_compress_on);

            if(lineItems.size() < 3){
                continue;
            }

            int type = boost::lexical_cast<int>(lineItems[1]);

            if(m_typeToCharge.size() < static_cast<size_t>(type + 1)){
                m_typeToCharge.resize(type + 1);
            }

            m_typeToCharge[type] = boost::lexical_cast<chemkit::Real>(lineItems[2]);
        }
    }

    return true;
}
