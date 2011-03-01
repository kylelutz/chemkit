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

#include "mcdlreader.h"

// === McdlReader ========================================================== //
// --- Construction and Destruction ---------------------------------------- //
McdlReader::McdlReader()
{
    p = 0;
    m_formula = 0;
    m_molecule = 0;
}

McdlReader::~McdlReader()
{
}

// --- Reading ------------------------------------------------------------- //
bool McdlReader::read(const QString &formula, chemkit::Molecule *molecule)
{
    return read(formula.toAscii().constData(), molecule);
}

bool McdlReader::read(const char *formula, chemkit::Molecule *molecule)
{
    bool ok;

    p = formula;
    m_formula = formula;
    m_molecule = molecule;

    // read composition module
    ok = readCompositionModule();
    if(!ok){
        return false;
    }

    // check for start of connectivity module
    if(*p == '['){
        p++;
    }
    else{
        return true;
    }

    // read connection module
    ok = readConnectionModule();
    if(!ok){
        return false;
    }

    return true;
}

// --- Error Handling ------------------------------------------------------ //
void McdlReader::setErrorString(const QString &error)
{
    m_errorString = error;
}

QString McdlReader::errorString() const
{
    return m_errorString;
}

// --- Internal Methods ---------------------------------------------------- //
bool McdlReader::readCompositionModule()
{
    chemkit::Atom *atom = 0;
    int quantity = 1;

    // read composition module
    for(;;){
        // check for end of string
        if(*p == '\0'){
            return true;
        }
        // check for '[' which signals start of connectivity module
        else if(*p == '['){
            addFragmentCopies(atom, quantity);

            break;
        }
        // check for ';' which signals start of next structural fragment
        else if(*p == ';'){
            addFragmentCopies(atom, quantity);

            atom = 0;
            quantity = 1;

            p++;
        }
        // check for quantity
        else if(isdigit(*p)){
            quantity = readNumber(&p);
        }
        // read element
        else if(isupper(*p)){
            if(atom){
                // add new terminal atom
                int atomicNumber = readElement(&p);

                chemkit::Atom *terminalAtom = m_molecule->addAtom(atomicNumber);
                if(!terminalAtom){
                    setErrorString(QString("Invalid terminal element in formula: %1").arg(atomicNumber));
                    return false;
                }

                m_molecule->addBond(atom, terminalAtom);
            }
            else{
                // add root atom
                int atomicNumber = readElement(&p);

                atom = m_molecule->addAtom(atomicNumber);
                if(!atom){
                    setErrorString(QString("Invalid element in formula: %1").arg(atomicNumber));
                    return false;
                }

                m_fragments.append(atom);
            }

            p++;
        }
        else{
            setErrorString(QString("Invalid character in formula: %1").arg(*p));
            return false;
        }
    }

    return true;
}

bool McdlReader::readConnectionModule()
{
    QList<int> connections;
    int fragment = 0;

    // read connectivity module
    for(;;){
        // check for end of string
        if(*p == '\0'){
            return true;
        }
        // check for end of connectivity module
        else if(*p == ']'){
            addFragmentConnections(connections, fragment);

            break;
        }
        else if(*p == ';'){
            addFragmentConnections(connections, fragment);

            connections.clear();
            fragment++;

            p++;
        }
        // read number
        else if(isdigit(*p)){
            connections.append(readNumber(&p));
        }
        else if(*p == ','){
            p++;
        }
        else{
            setErrorString(QString("Invalid character in formula: %1").arg(*p));
            return false;
        }
    }

    return true;
}

int McdlReader::readNumber(const char **p)
{
    int number = *(*p)++ - '0';

    while(isdigit(**p)){
        number = number * 10 + *(*p)++ - '0';
    }

    return number;
}

int McdlReader::readElement(const char **p)
{
    char first = **p;
    char second = *(*p+1);

    if(islower(second)){
        (*p)++;
        return chemkit::Element::atomicNumber(*p-1, 2);
    }
    else{
        return chemkit::Element::atomicNumber(first);
    }
}

// add 'quantity' number of copies of the fragment
void McdlReader::addFragmentCopies(chemkit::Atom *atom, int quantity)
{
    for(int i = 1; i < quantity; i++){
        chemkit::Atom *root = m_molecule->addAtomCopy(atom);
        m_fragments.append(root);

        foreach(chemkit::Atom *neighbor, atom->neighbors()){
            chemkit::Atom *terminalAtom  = m_molecule->addAtomCopy(neighbor);
            m_molecule->addBond(root, terminalAtom);
        }
    }
}

/// adds inter-fragment bonds for fragment
void McdlReader::addFragmentConnections(const QList<int> &connections, int fragment)
{
    chemkit::Atom *root = m_fragments.value(fragment);

    for(int i = 0; i < connections.size(); i++){
        chemkit::Atom *neighbor = m_fragments.value(connections.at(i) - 1);

        m_molecule->addBond(root, neighbor);
    }
}
