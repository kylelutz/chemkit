/******************************************************************************
**
** Copyright (C) 2009-2011 Kyle Lutz <kyle.r.lutz@gmail.com>
**
** This file is part of chemkit. For more information see
** <http://www.chemkit.org>.
**
** chemkit is free software: you can redistribute it and/or modify
** it under the terms of the GNU Lesser General Public License as published by
** the Free Software Foundation, either version 3 of the License, or
** (at your option) any later version.
**
** chemkit is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU Lesser General Public License for more details.
**
** You should have received a copy of the GNU Lesser General Public License
** along with chemkit. If not, see <http://www.gnu.org/licenses/>.
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
