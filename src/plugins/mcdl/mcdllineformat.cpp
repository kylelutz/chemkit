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

#include "mcdllineformat.h"

#include "mcdlreader.h"

// === McdlLineFormat ====================================================== //
// --- Construction and Destruction ---------------------------------------- //
McdlLineFormat::McdlLineFormat()
    : chemkit::LineFormat("mcdl")
{
}

McdlLineFormat::~McdlLineFormat()
{
}

// --- Input and Output ---------------------------------------------------- //
bool McdlLineFormat::read(const std::string &formula, chemkit::Molecule *molecule)
{
    McdlReader reader;

    bool ok = reader.read(formula, molecule);
    if(!ok){
        setErrorString(reader.errorString());
    }

    return ok;
}

std::string McdlLineFormat::write(const chemkit::Molecule *molecule)
{
    Q_UNUSED(molecule);

    setErrorString("MCDL write not supported.");

    return std::string();
}
