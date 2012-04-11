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

#include "mcdlreader.h"

#include <chemkit/atom.h>
#include <chemkit/foreach.h>

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
bool McdlReader::read(const std::string &formula, chemkit::Molecule *molecule)
{
    return read(formula.c_str(), molecule);
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
void McdlReader::setErrorString(const std::string &error)
{
    m_errorString = error;
}

std::string McdlReader::errorString() const
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
                chemkit::Element element = readElement(&p);

                chemkit::Atom *terminalAtom = m_molecule->addAtom(element);
                if(!terminalAtom){
                    setErrorString("Invalid terminal element in formula");
                    return false;
                }

                m_molecule->addBond(atom, terminalAtom);
            }
            else{
                // add root atom
                chemkit::Element element = readElement(&p);

                atom = m_molecule->addAtom(element);
                if(!atom){
                    setErrorString("Invalid element in formula");
                    return false;
                }

                m_fragments.push_back(atom);
            }

            p++;
        }
        else{
            setErrorString("Invalid character in formula");
            return false;
        }
    }

    return true;
}

bool McdlReader::readConnectionModule()
{
    std::vector<int> connections;
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
            connections.push_back(readNumber(&p));
        }
        else if(*p == ','){
            p++;
        }
        else{
            setErrorString("Invalid character in formula");
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

chemkit::Element McdlReader::readElement(const char **p)
{
    char first = **p;
    char second = *(*p+1);

    if(islower(second)){
        (*p)++;
        return chemkit::Element::fromSymbol(*p-1, 2);
    }
    else{
        return chemkit::Element::fromSymbol(first);
    }
}

// add 'quantity' number of copies of the fragment
void McdlReader::addFragmentCopies(chemkit::Atom *atom, int quantity)
{
    for(int i = 1; i < quantity; i++){
        chemkit::Atom *root = m_molecule->addAtomCopy(atom);
        m_fragments.push_back(root);

        foreach(chemkit::Atom *neighbor, atom->neighbors()){
            chemkit::Atom *terminalAtom  = m_molecule->addAtomCopy(neighbor);
            m_molecule->addBond(root, terminalAtom);
        }
    }
}

/// adds inter-fragment bonds for fragment
void McdlReader::addFragmentConnections(const std::vector<int> &connections, int fragment)
{
    chemkit::Atom *root = m_fragments[fragment];

    for(size_t i = 0; i < connections.size(); i++){
        chemkit::Atom *neighbor = m_fragments[connections[i] - 1];

        m_molecule->addBond(root, neighbor);
    }
}
