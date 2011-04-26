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

#include "mmffatom.h"
#include "mmffplugin.h"
#include "mmffforcefield.h"
#include "mmffparametersdata.h"

#include <chemkit/atom.h>
#include <chemkit/bond.h>
#include <chemkit/ring.h>
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
    d->deref();
}

// --- Parameters ---------------------------------------------------------- //
std::string MmffParameters::fileName() const
{
    return m_fileName;
}

bool MmffParameters::read(const std::string &fileName)
{
    // delete old parameters data
    if(d){
        d->deref();
        d = 0;
    }

    // try to load cached parameters
    MmffPlugin *mmffPlugin = static_cast<MmffPlugin *>(chemkit::PluginManager::instance()->plugin("mmff"));
    if(mmffPlugin){
        d = mmffPlugin->parameters(QString::fromStdString(fileName));

        if(d){
            d->ref();
            return true;
        }
    }

    // create new parameters data if we don't have a cached one
    if(!d){
        d = new MmffParametersData;
    }

    QFile file(QString::fromStdString(fileName));
    if(!file.open(QFile::ReadOnly)){
        setErrorString(file.errorString());
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

    while(!file.atEnd()){
        QString line = file.readLine();

        // lines that start with '$' indicate a new section
        if(line.startsWith("$")){
            section++;

            if(section == End){
                break;
            }
        }

        // lines starting with '#' are comments
        else if(line.startsWith("#")){
            continue;
        }

        // read data from line
        else{
            QStringList data = line.split(" ", QString::SkipEmptyParts);
            if(data.isEmpty() || data.size() < 2){
                continue;
            }

            if(section == BondStrech){
                int bondType = data.value(0).toInt();
                int typeA = data.value(1).toInt();
                int typeB = data.value(2).toInt();

                int index = calculateBondStrechIndex(bondType, typeA, typeB);

                MmffBondStrechParameters *parameters = new MmffBondStrechParameters;
                parameters->kb = data.value(3).toDouble();
                parameters->r0 = data.value(4).toDouble();

                d->bondStrechParameters[index] = parameters;
            }
            else if(section == EmpiricalBondStrech){
            }
            else if(section == AngleBend){
                int angleType = data.value(0).toInt();
                int typeA = data.value(1).toInt();
                int typeB = data.value(2).toInt();
                int typeC = data.value(3).toInt();

                int index = calculateAngleBendIndex(angleType, typeA, typeB, typeC);

                MmffAngleBendParameters *parameters = new MmffAngleBendParameters;
                parameters->ka = data.value(4).toDouble();
                parameters->theta0 = data.value(5).toDouble();

                d->angleBendParameters[index] = parameters;
            }
            else if(section == StrechBend){
                int strechBendType = data.value(0).toInt();
                int typeA = data.value(1).toInt();
                int typeB = data.value(2).toInt();
                int typeC = data.value(3).toInt();

                int index = calculateStrechBendIndex(strechBendType, typeA, typeB, typeC);

                MmffStrechBendParameters *parameters = new MmffStrechBendParameters;
                parameters->kba_ijk = data.value(4).toDouble();
                parameters->kba_kji = data.value(5).toDouble();

                d->strechBendParameters[index] = parameters;
            }
            else if(section == DefaultStrechBend){
                MmffDefaultStrechBendParameters *parameters = new MmffDefaultStrechBendParameters;
                parameters->rowA = data.value(0).toInt();
                parameters->rowB = data.value(1).toInt();
                parameters->rowC = data.value(2).toInt();
                parameters->parameters.kba_ijk = data.value(3).toDouble();
                parameters->parameters.kba_kji = data.value(4).toDouble();
                d->defaultStrechBendParameters.append(parameters);
            }
            else if(section == OutOfPlaneBending){
                int typeA = data.value(0).toInt();
                int typeB = data.value(1).toInt();
                int typeC = data.value(2).toInt();
                int typeD = data.value(3).toInt();

                int index = calculateOutOfPlaneBendingIndex(typeA, typeB, typeC, typeD);

                MmffOutOfPlaneBendingParameters *parameters = new MmffOutOfPlaneBendingParameters;
                parameters->koop = data.value(4).toDouble();

                d->outOfPlaneBendingParameters[index] = parameters;
            }
            else if(section == Torsion){
                int torsionType = data.value(0).toInt();
                int typeA = data.value(1).toInt();
                int typeB = data.value(2).toInt();
                int typeC = data.value(3).toInt();
                int typeD = data.value(4).toInt();

                int index = calculateTorsionIndex(torsionType, typeA, typeB, typeC, typeD);

                MmffTorsionParameters *parameters = new MmffTorsionParameters;
                parameters->V1 = data.value(5).toDouble();
                parameters->V2 = data.value(6).toDouble();
                parameters->V3 = data.value(7).toDouble();

                d->torsionParameters[index] = parameters;
            }
            else if(section == VanDerWaals){
                int type = data.value(0).toInt();
                if(type > MaxAtomType)
                    continue;

                MmffVanDerWaalsParameters *parameters = new MmffVanDerWaalsParameters;
                parameters->alpha = data.value(1).toDouble();
                parameters->N = data.value(2).toDouble();
                parameters->A = data.value(3).toDouble();
                parameters->G = data.value(4).toDouble();
                parameters->DA = data.value(5).at(0).toAscii();

                d->vanDerWaalsParameters[type] = parameters;
            }
            else if(section == Charge){
                MmffChargeParameters *parameters = new MmffChargeParameters;
                parameters->bondType = data.value(0).toInt();
                parameters->typeA = data.value(1).toInt();
                parameters->typeB = data.value(2).toInt();
                parameters->bci = data.value(3).toDouble();
                d->chargeParameters.append(parameters);
            }
            else if(section == PartialCharge){
                int type = data.value(1).toInt();
                if(type > MaxAtomType)
                    continue;

                MmffPartialChargeParameters *parameters = new MmffPartialChargeParameters;
                parameters->pbci = data.value(2).toDouble();
                parameters->fcadj = data.value(3).toDouble();

                d->partialChargeParameters[type] = parameters;
            }
        }
    }

    // store parameters in the cache
    if(mmffPlugin){
        mmffPlugin->storeParameters(QString::fromStdString(fileName), d);
    }

    return true;
}

const MmffBondStrechParameters* MmffParameters::bondStrechParameters(const MmffAtom *a, const MmffAtom *b) const
{
    int typeA = a->typeNumber();
    int typeB = b->typeNumber();
    int bondType = calculateBondType(a, b);

    return bondStrechParameters(bondType, typeA, typeB);
}

const MmffAngleBendParameters* MmffParameters::angleBendParameters(const MmffAtom *a, const MmffAtom *b, const MmffAtom *c) const
{
    int typeA = qMin(a->typeNumber(), c->typeNumber());
    int typeB = b->typeNumber();
    int typeC = qMax(a->typeNumber(), c->typeNumber());
    int angleType = calculateAngleType(a, b, c);

    return angleBendParameters(angleType, typeA, typeB, typeC);
}

const MmffStrechBendParameters* MmffParameters::strechBendParameters(const MmffAtom *a, const MmffAtom *b, const MmffAtom *c) const
{
    int typeA = a->typeNumber();
    int typeB = b->typeNumber();
    int typeC = c->typeNumber();
    int strechBendType = calculateStrechBendType(a, b, c);

    return strechBendParameters(strechBendType, typeA, typeB, typeC);
}

const MmffStrechBendParameters* MmffParameters::defaultStrechBendParameters(const MmffAtom *a, const MmffAtom *b, const MmffAtom *c) const
{
    return defaultStrechBendParameters(a->period()-1, b->period()-1, c->period()-1);
}

const MmffOutOfPlaneBendingParameters* MmffParameters::outOfPlaneBendingParameters(const MmffAtom *a, const MmffAtom *b, const MmffAtom *c, const MmffAtom *d) const
{
    int typeA = a->typeNumber();
    int typeB = b->typeNumber();
    int typeC = c->typeNumber();
    int typeD = d->typeNumber();

    const MmffOutOfPlaneBendingParameters *parameters = outOfPlaneBendingParameters(typeA, typeB, typeC, typeD);
    if(parameters)
        return parameters;

    // step down 3-2-3-3
    parameters = outOfPlaneBendingParameters(equivalentType(a, 3), typeB, equivalentType(c, 3), equivalentType(d, 3));
    if(parameters)
        return parameters;

    // step down 4-2-4-4
    parameters = outOfPlaneBendingParameters(equivalentType(a, 4), typeB, equivalentType(c, 4), equivalentType(d, 4));
    if(parameters)
        return parameters;

    // step down 5-2-5-5
    parameters = outOfPlaneBendingParameters(equivalentType(a, 5), typeB, equivalentType(c, 5), equivalentType(d, 5));
    if(parameters)
        return parameters;

    return 0;
}

const MmffTorsionParameters* MmffParameters::torsionParameters(const MmffAtom *a, const MmffAtom *b, const MmffAtom *c, const MmffAtom *d) const
{
    int typeA = a->typeNumber();
    int typeB = b->typeNumber();
    int typeC = c->typeNumber();
    int typeD = d->typeNumber();
    int torsionType = calculateTorsionType(a, b, c, d);

    const MmffTorsionParameters *parameters = torsionParameters(torsionType, typeA, typeB, typeC, typeD);
    if(parameters)
        return parameters;

    // step down 3-2-2-5
    parameters = torsionParameters(torsionType, equivalentType(a, 3), typeB, typeC, equivalentType(d, 5));
    if(parameters)
        return parameters;

    // step down 5-2-2-3
    parameters = torsionParameters(torsionType, equivalentType(a, 5), typeB, typeC, equivalentType(d, 3));
    if(parameters)
        return parameters;

    // step down 5-2-2-5
    parameters = torsionParameters(torsionType, equivalentType(a, 5), typeB, typeC, equivalentType(d, 5));
    if(parameters)
        return parameters;

    parameters = torsionParameters(0, equivalentType(a, 5), typeB, typeC, equivalentType(d, 5));
    if(parameters)
        return parameters;

    return 0;
}

const MmffVanDerWaalsParameters* MmffParameters::vanDerWaalsParameters(const MmffAtom *atom) const
{
    int type = atom->typeNumber();

    return d->vanDerWaalsParameters.value(type, 0);
}

const MmffAtomParameters* MmffParameters::atomParameters(int type) const
{
    if(type < 1 || type > 99)
        return 0;

    return &AtomParameters[type-1];
}

const MmffAtomParameters* MmffParameters::atomParameters(const MmffAtom *atom) const
{
    return atomParameters(atom->typeNumber());
}

const MmffChargeParameters* MmffParameters::chargeParameters(const chemkit::Atom *a, int typeA, const chemkit::Atom *b, int typeB) const
{
    int bondType = calculateBondType(a->bondTo(b), typeA, typeB);

    foreach(const MmffChargeParameters *parameters, d->chargeParameters){
        if(parameters->bondType == bondType &&
           parameters->typeA == typeA &&
           parameters->typeB == typeB){
            return parameters;
        }
    }

    return 0;
}

const MmffChargeParameters* MmffParameters::chargeParameters(const MmffAtom *a, const MmffAtom *b) const
{
    return chargeParameters(a->atom(), a->typeNumber(), b->atom(), b->typeNumber());
}

const MmffPartialChargeParameters* MmffParameters::partialChargeParameters(int type) const
{
    return d->partialChargeParameters.value(type, 0);
}

const MmffPartialChargeParameters* MmffParameters::partialChargeParameters(const MmffAtom *atom) const
{
    return partialChargeParameters(atom->typeNumber());
}

// --- Internal Methods ---------------------------------------------------- //
const MmffBondStrechParameters* MmffParameters::bondStrechParameters(int bondType, int typeA, int typeB) const
{
    if(typeA > typeB)
        qSwap(typeA, typeB);

    int index = calculateBondStrechIndex(bondType, typeA, typeB);

    return d->bondStrechParameters.value(index, 0);
}

const MmffBondStrechParameters* MmffParameters::empiricalBondStrechParameters(int atomicNumberA, int atomicNumberB) const
{
    Q_UNUSED(atomicNumberA);
    Q_UNUSED(atomicNumberB);

    return 0;
}

const MmffAngleBendParameters* MmffParameters::angleBendParameters(int angleType, int typeA, int typeB, int typeC) const
{
    if(typeA > typeC)
        qSwap(typeA, typeC);

    int index = calculateAngleBendIndex(angleType, typeA, typeB, typeC);

    return d->angleBendParameters.value(index, 0);
}

const MmffStrechBendParameters* MmffParameters::strechBendParameters(int strechBendType, int typeA, int typeB, int typeC) const
{
    int index = calculateStrechBendIndex(strechBendType, typeA, typeB, typeC);

    return d->strechBendParameters.value(index, 0);
}

const MmffStrechBendParameters* MmffParameters::defaultStrechBendParameters(int rowA, int rowB, int rowC) const
{
    foreach(const MmffDefaultStrechBendParameters *parameters, d->defaultStrechBendParameters){
        if(parameters->rowA == rowA &&
           parameters->rowB == rowB &&
           parameters->rowC == rowC){
            return &parameters->parameters;
        }
    }

    return 0;
}

const MmffOutOfPlaneBendingParameters* MmffParameters::outOfPlaneBendingParameters(int typeA, int typeB, int typeC, int typeD) const
{
    if(typeA > typeC)
        qSwap(typeA, typeD);

    if(typeC > typeD)
        qSwap(typeC, typeD);

    int index = calculateOutOfPlaneBendingIndex(typeA, typeB, typeC, typeD);

    return d->outOfPlaneBendingParameters.value(index, 0);
}

const MmffTorsionParameters* MmffParameters::torsionParameters(int torsionType, int typeA, int typeB, int typeC, int typeD) const
{
    if(typeB > typeC){
        qSwap(typeB, typeC);
        qSwap(typeA, typeD);
    }
    else if(typeB == typeC && typeA > typeD){
        qSwap(typeA, typeD);
    }

    int index = calculateTorsionIndex(torsionType, typeA, typeB, typeC, typeD);

    return d->torsionParameters.value(index, 0);
}

int MmffParameters::calculateBondType(const chemkit::Bond *bond, int typeA, int typeB) const
{
    const MmffAtomParameters *pa = atomParameters(typeA);
    const MmffAtomParameters *pb = atomParameters(typeB);
    if(!pa || !pb){
        return 0;
    }

    if(bond->order() == chemkit::Bond::Single && !MmffForceField::isAromatic(bond)){
        if(pa->sbmb && pb->sbmb){
            return 1;
        }
        else if(pa->arom && pb->arom){
            return 1;
        }
    }

    return 0;
}

int MmffParameters::calculateBondType(const MmffAtom *a, const MmffAtom *b) const
{
    const chemkit::Bond *bond = a->atom()->bondTo(b->atom());

    return calculateBondType(bond, a->typeNumber(), b->typeNumber());
}

int MmffParameters::calculateAngleType(const MmffAtom *a, const MmffAtom *b, const MmffAtom *c) const
{
    int bondTypeAB = calculateBondType(a, b);
    int bondTypeBC = calculateBondType(b, c);
    int bondTypeSum = bondTypeAB + bondTypeBC;

    bool inThreeMemberedRing = false;
    bool inFourMemberedRing = false;

    if(a->atom()->isBondedTo(c->atom())){
        inThreeMemberedRing = true;
    }
    else{
        foreach(const chemkit::Atom *neighbor, a->atom()->neighbors()){
            if(neighbor == b->atom())
                continue;

            if(neighbor->isBondedTo(c->atom()))
                inFourMemberedRing = true;
        }
    }

    if(inThreeMemberedRing){
        if(bondTypeSum == 1){
            return 5;
        }
        else if(bondTypeSum == 2){
            return 6;
        }
        else{
            return 3;
        }
    }
    else if(inFourMemberedRing){
        if(bondTypeSum == 1){
            return 7;
        }
        else if(bondTypeSum == 2){
            return 8;
        }
        else{
            return 4;
        }
    }
    else if(bondTypeSum == 1){
        return 1;
    }
    else if(bondTypeSum == 2){
        return 2;
    }
    else{
        return 0;
    }
}

int MmffParameters::calculateStrechBendType(const MmffAtom *a, const MmffAtom *b, const MmffAtom *c) const
{
    int bondTypeAB = calculateBondType(a, b);
    int bondTypeBC = calculateBondType(b, c);
    int angleType = calculateAngleType(a, b, c);

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

int MmffParameters::calculateTorsionType(const MmffAtom *a, const MmffAtom *b, const MmffAtom *c, const MmffAtom *d) const
{
    int bondTypeAB = calculateBondType(a, b);
    int bondTypeBC = calculateBondType(b, c);
    int bondTypeCD = calculateBondType(c, d);

    bool inFourMemberedRing = false;
    bool inFiveMemberedRing = false;

    if(a->atom()->isBondedTo(d->atom())){
        inFourMemberedRing = true;
    }

    foreach(const chemkit::Ring *ring, a->atom()->rings()){
        if(ring->size() == 5){
            if(ring->contains(b->atom()) && ring->contains(c->atom()) && ring->contains(d->atom())){
                if(!MmffForceField::isAromatic(ring)){
                    inFiveMemberedRing = true;
                    break;
                }
            }
        }
    }

    if(inFourMemberedRing){
        return 4;
    }
    else if(inFiveMemberedRing){
        return 5;
    }
    else if(bondTypeBC == 1){
        return 1;
    }
    else if(bondTypeAB == 1 || bondTypeCD == 1){
        return 2;
    }
    else{
        return 0;
    }
}

int MmffParameters::equivalentType(const MmffAtom *atom, int level) const
{
    if(level < 3){
        return atom->typeNumber();
    }

    for(int i = 0; i < EquivalentTypesCount; i++){
        if(EquivalentTypes[i][0] == atom->typeNumber()){
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
void MmffParameters::setErrorString(const QString &errorString)
{
    m_errorString = errorString;
}

QString MmffParameters::errorString() const
{
    return m_errorString;
}
