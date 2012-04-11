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

#include "spacedformulalineformat.h"

#include <map>

#include <chemkit/atom.h>
#include <chemkit/foreach.h>
#include <chemkit/molecule.h>

SpacedFormulaLineFormat::SpacedFormulaLineFormat()
    : chemkit::LineFormat("spaced-formula")
{
}

std::string SpacedFormulaLineFormat::write(const chemkit::Molecule *molecule)
{
    // a map of atomic symbols to their quantity
    std::map<std::string, size_t> composition;
    foreach(const chemkit::Atom *atom, molecule->atoms()){
        composition[atom->symbol()]++;
    }

    std::stringstream stream;

    if(composition.count("C") != 0){
        stream << "C ";
        stream << composition["C"] << " ";
        composition.erase("C");

        if(composition.count("H") != 0){
            stream << "H ";
            stream << composition["H"] << " ";
            composition.erase("H");
        }
    }

    std::map<std::string, size_t>::iterator iter;
    for(iter = composition.begin(); iter != composition.end(); ++iter){
        stream << iter->first << " ";
        stream << iter->second << " ";
    }

    std::string formula = stream.str();

    // remove trailing space
    if(!formula.empty()){
        formula.erase(formula.end() - 1);
    }

    return formula;
}
