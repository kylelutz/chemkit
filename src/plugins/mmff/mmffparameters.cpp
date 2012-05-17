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

#include "mmffparameters.h"

#include <fstream>

#include <boost/make_shared.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>

#include "mmffplugin.h"
#include "mmffatomtyper.h"
#include "mmffforcefield.h"
#include "mmffparametersdata.h"

#include <chemkit/atom.h>
#include <chemkit/bond.h>
#include <chemkit/foreach.h>
#include <chemkit/pluginmanager.h>

namespace {

// --- Atom Properties ----------------------------------------------------- //
const struct MmffAtomParameters AtomParameters[] = {
    {6, 4, 4, 0, 0, 0, 0, 0},
    {6, 3, 4, 0, 2, 0, 0, 1},
    {6, 3, 4, 0, 2, 0, 0, 1},
    {6, 2, 4, 0, 3, 0, 1, 1},
    {1, 1, 1, 0, 0, 0, 0, 0},
    {8, 2, 2, 1, 0, 0, 0, 0},
    {8, 1, 2, 0, 2, 0, 0, 0},
    {7, 3, 3, 1, 0, 0, 0, 0},
    {7, 2, 3, 0, 2, 0, 0, 1},
    {7, 3, 3, 1, 1, 0, 0, 0},
    {9, 1, 1, 1, 0, 0, 0, 0},
    {17, 1, 1, 1, 0, 0, 0, 0},
    {35, 1, 1, 1, 0, 0, 0, 0},
    {53, 1, 1, 1, 0, 0, 0, 0},
    {16, 2, 2, 1, 0, 0, 0, 0},
    {16, 1, 2, 0, 2, 0, 0, 0},
    {16, 3, 4, 0, 2, 0, 0, 0},
    {16, 4, 4, 0, 0, 0, 0, 0},
    {14, 4, 4, 0, 0, 0, 0, 0},
    {6, 4, 4, 0, 0, 0, 0, 0},
    {1, 1, 1, 0, 0, 0, 0, 0},
    {6, 4, 4, 0, 0, 0, 0, 0},
    {1, 1, 1, 0, 0, 0, 0, 0},
    {1, 1, 1, 0, 0, 0, 0, 0},
    {15, 4, 4, 0, 0, 0, 0, 0},
    {15, 3, 3, 1, 0, 0, 0, 0},
    {1, 1, 1, 0, 0, 0, 0, 0},
    {1, 1, 1, 0, 0, 0, 0, 0},
    {1, 1, 1, 0, 0, 0, 0, 0},
    {6, 3, 4, 0, 2, 0, 0, 1},
    {1, 1, 1, 0, 0, 0, 0, 0},
    {8, 1, 12, 1, 1, 0, 0, 0},
    {1, 1, 1, 0, 0, 0, 0, 0},
    {7, 4, 4, 0, 0, 0, 0, 0},
    {8, 1, 1, 1, 1, 0, 0, 0},
    {1, 1, 1, 0, 0, 0, 0, 0},
    {6, 3, 4, 0, 2, 1, 0, 1},
    {7, 2, 3, 0, 2, 1, 0, 0},
    {7, 3, 3, 1, 1, 1, 0, 1},
    {7, 3, 3, 1, 0, 0, 0, 0},
    {6, 3, 4, 0, 1, 0, 0, 0},
    {7, 1, 3, 0, 3, 0, 0, 0},
    {7, 3, 3, 1, 0, 0, 0, 0},
    {16, 2, 2, 1, 1, 1, 0, 0},
    {7, 3, 4, 0, 2, 0, 0, 0},
    {7, 2, 3, 0, 2, 0, 0, 0},
    {7, 1, 2, 0, 2, 0, 0, 0},
    {7, 2, 2, 0, 0, 0, 0, 0},
    {8, 3, 3, 0, 0, 0, 0, 0},
    {1, 1, 1, 0, 0, 0, 0, 0},
    {8, 2, 3, 0, 2, 0, 0, 0},
    {1, 1, 1, 0, 0, 0, 0, 0},
    {7, 2, 4, 0, 2, 0, 1, 0},
    {7, 3, 4, 0, 2, 0, 0, 1},
    {7, 3, 34, 0, 1, 0, 0, 0},
    {7, 3, 34, 0, 1, 0, 0, 0},
    {6, 3, 4, 0, 2, 0, 0, 1},
    {7, 3, 4, 0, 1, 1, 0, 1},
    {8, 2, 2, 1, 1, 1, 0, 0},
    {6, 1, 3, 0, 3, 0, 0, 0},
    {7, 2, 4, 0, 3, 0, 1, 0},
    {7, 2, 2, 1, 0, 0, 0, 0},
    {6, 3, 4, 0, 2, 1, 0, 1},
    {6, 3, 4, 0, 2, 1, 0, 1},
    {7, 2, 3, 0, 2, 1, 0, 0},
    {7, 2, 3, 0, 2, 1, 0, 0},
    {7, 3, 4, 0, 2, 0, 0, 1},
    {7, 4, 4, 0, 0, 0, 0, 0},
    {7, 3, 4, 0, 1, 1, 0, 0},
    {8, 2, 2, 1, 0, 0, 0, 0},
    {1, 1, 1, 0, 0, 0, 0, 0},
    {16, 1, 1, 1, 1, 0, 0, 0},
    {16, 3, 3, 0, 0, 0, 0, 0},
    {16, 2, 4, 0, 2, 0, 0, 0},
    {15, 2, 3, 0, 2, 0, 0, 1},
    {7, 2, 2, 1, 0, 0, 0, 0},
    {17, 4, 4, 0, 0, 0, 0, 0},
    {6, 3, 4, 0, 2, 1, 0, 1},
    {7, 2, 3, 0, 2, 1, 0, 0},
    {6, 3, 4, 0, 2, 0, 0, 1},
    {7, 3, 4, 0, 1, 1, 0, 1},
    {7, 3, 4, 0, 1, 1, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 0},
    {26, 0, 0, 0, 0, 0, 0, 0},
    {26, 0, 0, 0, 0, 0, 0, 0},
    {9, 0, 0, 0, 0, 0, 0, 0},
    {17, 0, 0, 0, 0, 0, 0, 0},
    {35, 0, 0, 0, 0, 0, 0, 0},
    {3, 0, 0, 0, 0, 0, 0, 0},
    {11, 0, 0, 0, 0, 0, 0, 0},
    {19, 0, 0, 0, 0, 0, 0, 0},
    {30, 0, 0, 0, 0, 0, 0, 0},
    {20, 0, 0, 0, 0, 0, 0, 0},
    {29, 0, 0, 0, 0, 0, 0, 0},
    {29, 0, 0, 0, 0, 0, 0, 0},
    {12, 0, 0, 0, 0, 0, 0, 0},
};

// --- Equivalent Types ---------------------------------------------------- //
const int EquivalentTypes[][5] = {
    {1, 1, 1, 1, 0},
    {2, 2, 2, 1, 0},
    {3, 3, 3, 1, 0},
    {4, 4, 4, 1, 0},
    {5, 5, 5, 5, 0},
    {6, 6, 6, 6, 0},
    {7, 7, 7, 6, 0},
    {8, 8, 8, 8, 0},
    {9, 9, 9, 8, 0},
    {10, 10, 10, 8, 0},
    {11, 11, 11, 11, 0},
    {12, 12, 12, 12, 0},
    {13, 13, 13, 13, 0},
    {14, 14, 14, 14, 0},
    {15, 15, 15, 15, 0},
    {16, 16, 16, 15, 0},
    {17, 17, 17, 15, 0},
    {18, 18, 18, 15, 0},
    {19, 19, 19, 19, 0},
    {20, 20, 1, 1, 0},
    {21, 21, 21, 5, 0},
    {22, 22, 22, 1, 0},
    {23, 23, 23, 5, 0},
    {24, 24, 24, 5, 0},
    {25, 25, 25, 25, 0},
    {26, 26, 26, 25, 0},
    {27, 27, 28, 5, 0},
    {28, 28, 28, 5, 0},
    {29, 29, 29, 5, 0},
    {30, 30, 2, 1, 0},
    {31, 31, 31, 31, 0},
    {32, 32, 7, 6, 0},
    {33, 33, 21, 5, 0},
    {34, 34, 8, 8, 0},
    {35, 35, 6, 6, 0},
    {36, 36, 36, 5, 0},
    {37, 37, 2, 1, 0},
    {38, 38, 9, 8, 0},
    {39, 39, 10, 8, 0},
    {40, 40, 10, 8, 0},
    {41, 41, 3, 1, 0},
    {42, 42, 42, 8, 0},
    {43, 43, 10, 8, 0},
    {44, 44, 16, 15, 0},
    {45, 45, 10, 8, 0},
    {46, 46, 9, 8, 0},
    {47, 47, 42, 8, 0},
    {48, 48, 9, 8, 0},
    {49, 49, 6, 6, 0},
    {50, 50, 21, 5, 0},
    {51, 51, 7, 6, 0},
    {52, 52, 21, 5, 0},
    {53, 53, 42, 8, 0},
    {54, 54, 9, 8, 0},
    {55, 55, 10, 8, 0},
    {56, 56, 10, 8, 0},
    {57, 57, 2, 1, 0},
    {58, 58, 10, 8, 0},
    {59, 59, 6, 6, 0},
    {60, 60, 4, 1, 0},
    {61, 61, 42, 8, 0},
    {62, 62, 10, 8, 0},
    {63, 63, 2, 1, 0},
    {64, 64, 2, 1, 0},
    {65, 65, 9, 8, 0},
    {66, 66, 9, 8, 0},
    {67, 67, 9, 8, 0},
    {68, 68, 8, 8, 0},
    {69, 69, 9, 8, 0},
    {70, 70, 70, 70, 70},
    {71, 71, 5, 5, 0},
    {72, 72, 16, 15, 0},
    {73, 73, 18, 15, 0},
    {74, 74, 17, 15, 0},
    {75, 75, 26, 25, 0},
    {76, 76, 9, 8, 0},
    {77, 77, 12, 12, 0},
    {78, 78, 2, 1, 0},
    {79, 79, 9, 8, 0},
    {80, 80, 2, 1, 0},
    {81, 81, 10, 8, 0},
    {82, 82, 9, 8, 0},
    {87, 87, 87, 87, 87},
    {88, 88, 88, 88, 88},
    {89, 89, 89, 89, 89},
    {90, 90, 90, 90, 90},
    {91, 91, 91, 91, 91},
    {92, 92, 92, 92, 92},
    {93, 93, 93, 93, 93},
    {94, 94, 94, 94, 94},
    {95, 95, 95, 95, 95},
    {96, 96, 96, 96, 96},
    {97, 97, 97, 97, 97},
    {98, 98, 98, 98, 98},
    {99, 99, 99, 99, 99},
};

int EquivalentTypesCount = sizeof(EquivalentTypes) / sizeof(*EquivalentTypes);

} // end anonymous namespace

// --- Construction and Destruction ---------------------------------------- //
MmffParameters::MmffParameters()
    : d(new MmffParametersData)
{
}

MmffParameters::~MmffParameters()
{
}

// --- Parameters ---------------------------------------------------------- //
std::string MmffParameters::fileName() const
{
    return m_fileName;
}

bool MmffParameters::read(const std::string &fileName)
{
    // try to load cached parameters
    MmffPlugin *mmffPlugin = static_cast<MmffPlugin *>(chemkit::PluginManager::instance()->plugin("mmff"));
    if(mmffPlugin){
        d = mmffPlugin->parameters(fileName);

        if(d){
            return true;
        }
    }

    // create new parameters data if we don't have a cached one
    if(!d){
        d = boost::make_shared<MmffParametersData>();
    }

    std::ifstream file(fileName.c_str());
    if(!file.is_open()){
        setErrorString("Failed to open parameters file.");
        return false;
    }

    m_fileName = fileName;

    // section in file
    enum Section {
        BondStrech,
        EmpiricalBondStrech,
        AngleBend,
        StrechBend,
        DefaultStrechBend,
        OutOfPlaneBending,
        Torsion,
        VanDerWaals,
        Charge,
        PartialCharge,
        End
    };

    // first section is bond strech parameters
    int section = BondStrech;

    while(!file.eof()){
        std::string line;
        std::getline(file, line);
        boost::trim_left(line);

        // lines that start with '$' indicate a new section
        if(boost::starts_with(line, "$")){
            section++;

            if(section == End){
                break;
            }
        }

        // lines starting with '#' are comments
        else if(boost::starts_with(line, "#")){
            continue;
        }

        // read data from line
        else{
            std::vector<std::string> data;
            boost::split(data, line, boost::is_any_of(" "), boost::token_compress_on);
            if(data.empty() || data.size() < 2){
                continue;
            }

            if(section == BondStrech){
                int bondType = boost::lexical_cast<int>(data[0]);
                int typeA = boost::lexical_cast<int>(data[1]);
                int typeB = boost::lexical_cast<int>(data[2]);

                int index = calculateBondStrechIndex(bondType, typeA, typeB);

                MmffBondStrechParameters parameters;
                parameters.kb = boost::lexical_cast<chemkit::Real>(data[3]);
                parameters.r0 = boost::lexical_cast<chemkit::Real>(data[4]);
                d->bondStrechParameters[index] = parameters;
            }
            else if(section == EmpiricalBondStrech){
            }
            else if(section == AngleBend){
                int angleType = boost::lexical_cast<int>(data[0]);
                int typeA = boost::lexical_cast<int>(data[1]);
                int typeB = boost::lexical_cast<int>(data[2]);
                int typeC = boost::lexical_cast<int>(data[3]);

                int index = calculateAngleBendIndex(angleType, typeA, typeB, typeC);

                MmffAngleBendParameters parameters;
                parameters.ka = boost::lexical_cast<chemkit::Real>(data[4]);
                parameters.theta0 = boost::lexical_cast<chemkit::Real>(data[5]);
                d->angleBendParameters[index] = parameters;
            }
            else if(section == StrechBend){
                int strechBendType = boost::lexical_cast<int>(data[0]);
                int typeA = boost::lexical_cast<int>(data[1]);
                int typeB = boost::lexical_cast<int>(data[2]);
                int typeC = boost::lexical_cast<int>(data[3]);

                int index = calculateStrechBendIndex(strechBendType, typeA, typeB, typeC);

                MmffStrechBendParameters parameters;
                parameters.kba_ijk = boost::lexical_cast<chemkit::Real>(data[4]);
                parameters.kba_kji = boost::lexical_cast<chemkit::Real>(data[5]);
                d->strechBendParameters[index] = parameters;
            }
            else if(section == DefaultStrechBend){
                MmffDefaultStrechBendParameters parameters;
                parameters.rowA = boost::lexical_cast<int>(data[0]);
                parameters.rowB = boost::lexical_cast<int>(data[1]);
                parameters.rowC = boost::lexical_cast<int>(data[2]);
                parameters.parameters.kba_ijk = boost::lexical_cast<chemkit::Real>(data[3]);
                parameters.parameters.kba_kji = boost::lexical_cast<chemkit::Real>(data[4]);
                d->defaultStrechBendParameters.push_back(parameters);
            }
            else if(section == OutOfPlaneBending){
                int typeA = boost::lexical_cast<int>(data[0]);
                int typeB = boost::lexical_cast<int>(data[1]);
                int typeC = boost::lexical_cast<int>(data[2]);
                int typeD = boost::lexical_cast<int>(data[3]);

                int index = calculateOutOfPlaneBendingIndex(typeA, typeB, typeC, typeD);

                MmffOutOfPlaneBendingParameters parameters;
                parameters.koop = boost::lexical_cast<chemkit::Real>(data[4]);
                d->outOfPlaneBendingParameters[index] = parameters;
            }
            else if(section == Torsion){
                int torsionType = boost::lexical_cast<int>(data[0]);
                int typeA = boost::lexical_cast<int>(data[1]);
                int typeB = boost::lexical_cast<int>(data[2]);
                int typeC = boost::lexical_cast<int>(data[3]);
                int typeD = boost::lexical_cast<int>(data[4]);

                int index = calculateTorsionIndex(torsionType, typeA, typeB, typeC, typeD);

                MmffTorsionParameters parameters;
                parameters.V1 = boost::lexical_cast<chemkit::Real>(data[5]);
                parameters.V2 = boost::lexical_cast<chemkit::Real>(data[6]);
                parameters.V3 = boost::lexical_cast<chemkit::Real>(data[7]);
                d->torsionParameters[index] = parameters;
            }
            else if(section == VanDerWaals){
                int type = boost::lexical_cast<int>(data[0]);
                if(type > MaxAtomType)
                    continue;

                MmffVanDerWaalsParameters parameters;
                parameters.alpha = boost::lexical_cast<chemkit::Real>(data[1]);
                parameters.N = boost::lexical_cast<chemkit::Real>(data[2]);
                parameters.A = boost::lexical_cast<chemkit::Real>(data[3]);
                parameters.G = boost::lexical_cast<chemkit::Real>(data[4]);
                parameters.DA = data[5][0];
                d->vanDerWaalsParameters[type] = parameters;
            }
            else if(section == Charge){
                MmffChargeParameters parameters;
                parameters.bondType = boost::lexical_cast<int>(data[0]);
                parameters.typeA = boost::lexical_cast<int>(data[1]);
                parameters.typeB = boost::lexical_cast<int>(data[2]);
                parameters.bci = boost::lexical_cast<chemkit::Real>(data[3]);
                d->chargeParameters.push_back(parameters);
            }
            else if(section == PartialCharge){
                int type = boost::lexical_cast<int>(data[1]);
                if(type > MaxAtomType)
                    continue;

                MmffPartialChargeParameters parameters;
                parameters.pbci = boost::lexical_cast<chemkit::Real>(data[2]);
                parameters.fcadj = boost::lexical_cast<chemkit::Real>(data[3]);
                d->partialChargeParameters[type] = parameters;
            }
        }
    }

    // store parameters in the cache
    if(mmffPlugin){
        mmffPlugin->storeParameters(fileName, d);
    }

    return true;
}

const MmffVanDerWaalsParameters* MmffParameters::vanDerWaalsParameters(int type) const
{
    return &d->vanDerWaalsParameters[type];
}

const MmffAtomParameters* MmffParameters::atomParameters(int type) const
{
    if(type < 1 || type > 99)
        return 0;

    return &AtomParameters[type-1];
}

const MmffChargeParameters* MmffParameters::chargeParameters(const chemkit::Atom *a, int typeA, const chemkit::Atom *b, int typeB) const
{
    int bondType = calculateBondType(a->bondTo(b), typeA, typeB);

    foreach(const MmffChargeParameters &parameters, d->chargeParameters){
        if(parameters.bondType == bondType &&
           parameters.typeA == typeA &&
           parameters.typeB == typeB){
            return &parameters;
        }
    }

    return 0;
}

const MmffPartialChargeParameters* MmffParameters::partialChargeParameters(int type) const
{
    return &d->partialChargeParameters[type];
}

const MmffBondStrechParameters* MmffParameters::bondStrechParameters(int bondType, int typeA, int typeB) const
{
    if(typeA > typeB)
        std::swap(typeA, typeB);

    int index = calculateBondStrechIndex(bondType, typeA, typeB);

    std::map<int, MmffBondStrechParameters>::const_iterator iter = d->bondStrechParameters.find(index);
    if(iter == d->bondStrechParameters.end()){
        return 0;
    }

    return &iter->second;
}

const MmffBondStrechParameters* MmffParameters::empiricalBondStrechParameters(int atomicNumberA, int atomicNumberB) const
{
    CHEMKIT_UNUSED(atomicNumberA);
    CHEMKIT_UNUSED(atomicNumberB);

    return 0;
}

const MmffAngleBendParameters* MmffParameters::angleBendParameters(int angleType, int typeA, int typeB, int typeC) const
{
    if(typeA > typeC)
        std::swap(typeA, typeC);

    int index = calculateAngleBendIndex(angleType, typeA, typeB, typeC);

    std::map<int, MmffAngleBendParameters>::const_iterator iter = d->angleBendParameters.find(index);
    if(iter == d->angleBendParameters.end()){
        return 0;
    }

    return &iter->second;
}

const MmffStrechBendParameters* MmffParameters::strechBendParameters(int strechBendType, int typeA, int typeB, int typeC) const
{
    int index = calculateStrechBendIndex(strechBendType, typeA, typeB, typeC);

    std::map<int, MmffStrechBendParameters>::const_iterator iter = d->strechBendParameters.find(index);
    if(iter == d->strechBendParameters.end()){
        return 0;
    }

    return &iter->second;
}

const MmffStrechBendParameters* MmffParameters::defaultStrechBendParameters(int typeA, int typeB, int typeC) const
{
    int rowA = MmffAtomTyper::typeToElement(typeA).period() - 1;
    int rowB = MmffAtomTyper::typeToElement(typeB).period() - 1;
    int rowC = MmffAtomTyper::typeToElement(typeC).period() - 1;

    foreach(const MmffDefaultStrechBendParameters &parameters, d->defaultStrechBendParameters){
        if(parameters.rowA == rowA &&
           parameters.rowB == rowB &&
           parameters.rowC == rowC){
            return &parameters.parameters;
        }
    }

    return 0;
}

const MmffOutOfPlaneBendingParameters* MmffParameters::outOfPlaneBendingParameters(int typeA, int typeB, int typeC, int typeD) const
{
    if(typeA > typeC)
        std::swap(typeA, typeC);
    if(typeA > typeD)
        std::swap(typeA, typeD);
    if(typeC > typeD)
        std::swap(typeC, typeD);

    int index = calculateOutOfPlaneBendingIndex(typeA, typeB, typeC, typeD);

    std::map<int, MmffOutOfPlaneBendingParameters>::const_iterator iter;
    iter = d->outOfPlaneBendingParameters.find(index);
    if(iter != d->outOfPlaneBendingParameters.end()){
        return &iter->second;
    }

    // step down 3-2-3-3
    index = calculateOutOfPlaneBendingIndex(equivalentType(typeA, 3), typeB, equivalentType(typeC, 3), equivalentType(typeD, 3));
    iter = d->outOfPlaneBendingParameters.find(index);
    if(iter != d->outOfPlaneBendingParameters.end()){
        return &iter->second;
    }

    // step down 4-2-4-4
    index = calculateOutOfPlaneBendingIndex(equivalentType(typeA, 4), typeB, equivalentType(typeC, 4), equivalentType(typeD, 4));
    iter = d->outOfPlaneBendingParameters.find(index);
    if(iter != d->outOfPlaneBendingParameters.end()){
        return &iter->second;
    }

    // step down 5-2-5-5
    index = calculateOutOfPlaneBendingIndex(equivalentType(typeA, 5), typeB, equivalentType(typeC, 5), equivalentType(typeD, 5));
    iter = d->outOfPlaneBendingParameters.find(index);
    if(iter != d->outOfPlaneBendingParameters.end()){
        return &iter->second;
    }

    return 0;
}

const MmffTorsionParameters* MmffParameters::torsionParameters(int torsionType, int typeA, int typeB, int typeC, int typeD) const
{
    if(typeB > typeC){
        std::swap(typeB, typeC);
        std::swap(typeA, typeD);
    }
    else if(typeB == typeC && typeA > typeD){
        std::swap(typeA, typeD);
    }

    int index = calculateTorsionIndex(torsionType, typeA, typeB, typeC, typeD);

    std::map<int, MmffTorsionParameters>::const_iterator iter;
    iter = d->torsionParameters.find(index);
    if(iter != d->torsionParameters.end()){
        return &iter->second;
    }

    // step down 3-2-2-5
    index = calculateTorsionIndex(torsionType, equivalentType(typeA, 3), typeB, typeC, equivalentType(typeD, 5));
    iter = d->torsionParameters.find(index);
    if(iter != d->torsionParameters.end()){
        return &iter->second;
    }

    // step down 5-2-2-3
    index = calculateTorsionIndex(torsionType, equivalentType(typeA, 5), typeB, typeC, equivalentType(typeD, 3));
    iter = d->torsionParameters.find(index);
    if(iter != d->torsionParameters.end()){
        return &iter->second;
    }

    // step down 5-2-2-5
    index = calculateTorsionIndex(torsionType, equivalentType(typeA, 5), typeB, typeC, equivalentType(typeD, 5));
    iter = d->torsionParameters.find(index);
    if(iter != d->torsionParameters.end()){
        return &iter->second;
    }

    // step down 5-2-2-5 (with no torsion type)
    index = calculateTorsionIndex(0, equivalentType(typeA, 5), typeB, typeC, equivalentType(typeD, 5));
    iter = d->torsionParameters.find(index);
    if(iter != d->torsionParameters.end()){
        return &iter->second;
    }

    return 0;
}

// --- Static Methods ------------------------------------------------------ //
int MmffParameters::calculateBondType(const chemkit::Bond *bond, int typeA, int typeB)
{
    MmffParameters parameters;
    const MmffAtomParameters *parametersA = parameters.atomParameters(typeA);
    const MmffAtomParameters *parametersB = parameters.atomParameters(typeB);

    if(!parametersA || !parametersB){
        return 0;
    }

    if(bond->order() == chemkit::Bond::Single && !MmffAromaticityModel().isAromatic(bond)){
        if(parametersA->sbmb && parametersB->sbmb){
            return 1;
        }
        else if(parametersA->arom && parametersB->arom){
            return 1;
        }
    }

    return 0;
}

int MmffParameters::calculateStrechBendType(int bondTypeAB, int bondTypeBC, int angleType)
{
    if(angleType == 0){
        return 0;
    }
    else if(angleType == 1){
        if(bondTypeAB == 1){
            return 1;
        }
        else if(bondTypeBC == 1){
            return 2;
        }
        else{
            return 0;
        }
    }
    else if(angleType == 2){
        return 3;
    }
    else if(angleType == 3){
        return 5;
    }
    else if(angleType == 4){
        return 4;
    }
    else if(angleType == 5){
        if(bondTypeAB == 1){
            return 6;
        }
        else if(bondTypeBC == 1){
            return 7;
        }
        else{
            return 0;
        }
    }
    else if(angleType == 6){
        return 8;
    }
    else if(angleType == 7){
        if(bondTypeAB == 1){
            return 9;
        }
        else if(bondTypeBC == 1){
            return 10;
        }
        else{
            return 0;
        }
    }
    else if(angleType == 8){
        return 11;
    }
    else{
        return 0;
    }
}

// --- Internal Methods ---------------------------------------------------- //
int MmffParameters::equivalentType(int type, int level) const
{
    if(level < 3){
        return type;
    }

    for(int i = 0; i < EquivalentTypesCount; i++){
        if(EquivalentTypes[i][0] == type){
            return EquivalentTypes[i][level-1];
        }
    }

    return 0;
}

int MmffParameters::calculateBondStrechIndex(int bondType, int typeA, int typeB) const
{
    return 2 * (typeA * 136 + typeB) + bondType;
}

int MmffParameters::calculateAngleBendIndex(int angleType, int typeA, int typeB, int typeC) const
{
    return 9 * (typeB * (136*136) + typeA * 136 + typeC) + angleType;
}

int MmffParameters::calculateStrechBendIndex(int strechBendType, int typeA, int typeB, int typeC) const
{
    return 12 * (typeB * (136*136) + typeA * 136 + typeC) + strechBendType;
}

int MmffParameters::calculateOutOfPlaneBendingIndex(int typeA, int typeB, int typeC, int typeD) const
{
    return typeB * (136*136*136) + typeA * (136*136) + typeC * 136 + typeD;
}

int MmffParameters::calculateTorsionIndex(int torsionType, int typeA, int typeB, int typeC, int typeD) const
{
    return 6 * (typeB * (136*136*136) + typeC * (136*136) + typeA * 136 + typeD) + torsionType;
}

// --- Error Handling ------------------------------------------------------ //
void MmffParameters::setErrorString(const std::string &errorString)
{
    m_errorString = errorString;
}

std::string MmffParameters::errorString() const
{
    return m_errorString;
}
