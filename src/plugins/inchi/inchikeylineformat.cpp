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

#include "inchikeylineformat.h"

#include "../../3rdparty/inchi/inchi_api.h"

InchiKeyLineFormat::InchiKeyLineFormat()
    : chemkit::LineFormat("inchikey")
{
}

std::string InchiKeyLineFormat::write(const chemkit::Molecule *molecule)
{
    LineFormat *inchiLineFormat = chemkit::LineFormat::create("inchi");
    if(!inchiLineFormat){
        setErrorString("Failed to load the InChI line format.");
        return std::string();
    }

    std::string inchi = inchiLineFormat->write(molecule);
    if(inchi.empty()){
        setErrorString(inchiLineFormat->errorString());
        delete inchiLineFormat;
        return std::string();
    }

    delete inchiLineFormat;

    char inchiKey[28]; // 27 characters + null terminator

    int ret = GetStdINCHIKeyFromStdINCHI(inchi.c_str(), inchiKey);
    if(ret != INCHIKEY_OK){
        setErrorString(QString("InChIKey failed: the generator returned '%1'.").arg(ret).toStdString());
        return std::string();
    }

    return std::string(inchiKey);
}
