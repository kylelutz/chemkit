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

#include <QtCore>

// --- Construction and Destruction ---------------------------------------- //
OplsParameters::OplsParameters(const QString &fileName)
    : m_fileName(fileName)
{
    read(fileName);
}

OplsParameters::~OplsParameters()
{
}

// --- Properties ---------------------------------------------------------- //
void OplsParameters::setFileName(const QString &fileName)
{
    m_fileName = fileName;
}

QString OplsParameters::fileName() const
{
    return m_fileName;
}

// --- Parameters ---------------------------------------------------------- //
int OplsParameters::atomClass(int type) const
{
    return m_typeToClass.value(type);
}

QString OplsParameters::atomName(int type) const
{
    return m_typeToName.value(type);
}

chemkit::Float OplsParameters::partialCharge(int type) const
{
    return m_typeToCharge.value(type);
}

const OplsBondStrechParameters* OplsParameters::bondStrechParameters(int a, int b) const
{
    a = atomClass(a);
    b = atomClass(b);

    if(a > b){
        qSwap(a, b);
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
        qSwap(a, c);
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
        qSwap(b, c);
    }

    if(a > d){
        qSwap(a, d);
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
bool OplsParameters::read(const QString &fileName)
{
    QFile file(fileName);
    bool ok = file.open(QFile::ReadOnly);
    if(!ok){
        return false;
    }

    while(!file.atEnd()){
        QByteArray line = file.readLine();

        // atom parameters
        if(line.startsWith("atom")){
            QStringList lineItems = QString(line).split(' ', QString::SkipEmptyParts);

            if(lineItems.size() < 4){
                continue;
            }

            int type = lineItems[1].toInt();
            int klass = lineItems[2].toInt();
            QString name = lineItems[3];

            if(m_typeToClass.size() < type + 1)
                m_typeToClass.resize(type + 1);
            m_typeToClass[type] = klass;

            if(m_typeToName.size() < type + 1)
                m_typeToName.resize(type + 1);
            m_typeToName[type] = name;
        }
        // bond parameters
        else if(line.startsWith("bond")){
            QStringList lineItems = QString(line).split(' ', QString::SkipEmptyParts);

            if(lineItems.size() < 5){
                continue;
            }

            OplsBondStrechParameters p;
            p.typeA = lineItems[1].toInt();
            p.typeB = lineItems[2].toInt();
            p.kb = lineItems[3].toDouble();
            p.r0 = lineItems[4].toDouble();

            m_bondStrechParameters.append(p);
        }
        // angle parameters
        else if(line.startsWith("angle")){
            QStringList lineItems = QString(line).split(' ', QString::SkipEmptyParts);

            if(lineItems.size() < 6){
                continue;
            }

            OplsAngleBendParameters p;
            p.typeA = lineItems[1].toInt();
            p.typeB = lineItems[2].toInt();
            p.typeC = lineItems[3].toInt();
            p.ka = lineItems[4].toDouble();
            p.theta0 = lineItems[5].toDouble();

            m_angleBendParameters.append(p);
        }
        // torsion parameters
        else if(line.startsWith("torsion")){
            QStringList lineItems = QString(line).split(' ', QString::SkipEmptyParts);

            if(lineItems.size() < 14){
                continue;
            }

            OplsTorsionParameters p;
            p.typeA = lineItems[1].toInt();
            p.typeB = lineItems[2].toInt();
            p.typeC = lineItems[3].toInt();
            p.typeD = lineItems[4].toInt();
            p.v1 = lineItems[5].toDouble();
            p.v2 = lineItems[8].toDouble();
            p.v3 = lineItems[11].toDouble();

            m_torsionParameters.append(p);
        }
        // van der waals parameters
        else if(line.startsWith("vdw")){
            QStringList lineItems = QString(line).split(' ', QString::SkipEmptyParts);

            if(lineItems.size() < 4){
                continue;
            }

            int type = lineItems[1].toInt();
            if(m_vanDerWaalsParameters.size() < type + 1)
                m_vanDerWaalsParameters.resize(type + 1);

            OplsVanDerWaalsParameters p;
            p.sigma = lineItems[2].toDouble();
            p.epsilon = lineItems[3].toDouble();

            m_vanDerWaalsParameters[type] = p;
        }
        else if(line.startsWith("charge")){
            QStringList lineItems = QString(line).split(' ', QString::SkipEmptyParts);

            if(lineItems.size() < 3){
                continue;
            }

            int type = lineItems[1].toInt();

            if(m_typeToCharge.size() < type + 1){
                m_typeToCharge.resize(type + 1);
            }

            m_typeToCharge[type] = lineItems[2].toDouble();
        }
    }

    return true;
}
