/******************************************************************************
**
** Copyright (C) 2012 Kitware, Inc.
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

#include "chemjsonfileformat.h"

#include <boost/range/algorithm.hpp>

#include <chemkit/atom.h>
#include <chemkit/bond.h>
#include <chemkit/foreach.h>
#include <chemkit/molecule.h>
#include <chemkit/moleculefile.h>

#include "../../3rdparty/jsoncpp/json/json.h"

namespace {

// Returns a sanitized copy of the string which can then be
// safely inserted into a JSON document.
std::string sanitizeJsonString(std::string string)
{
    // remove any double-quotes
    string.erase(boost::remove(string, '"'), string.end());

    return string;
}

} // end anonymous namespace

// The ChemJsonFileFormat class implements reading and writing
// for files in the Chemical JSON file format.
//
// Specification: http://wiki.openchemistry.org/Chemical_JSON
ChemJsonFileFormat::ChemJsonFileFormat()
    : chemkit::MoleculeFileFormat("cjson")
{
}

bool ChemJsonFileFormat::read(std::istream &input, chemkit::MoleculeFile *file)
{
    // read file
    Json::Value root;
    Json::Reader reader;

    bool ok = reader.parse(input, root);
    if(!ok){
        setErrorString(reader.getFormattedErrorMessages());
        return false;
    }

    // check file type
    Json::Value version = root["chemical json"];
    if(version.empty()){
        setErrorString("Not a valid Chemical JSON file");
        return false;
    }

    // create molecule
    boost::shared_ptr<chemkit::Molecule> molecule(new chemkit::Molecule);

    // read name
    Json::Value name = root["name"];
    if(!name.empty() && name.type() == Json::stringValue){
        molecule->setName(name.asString());
    }

    // read elements
    Json::Value elements = root["atoms"]["elements"];
    for(Json::ArrayIndex i = 0; i < elements.size(); i++){
        molecule->addAtom(elements.get(i, 0).asInt());
    }

    // read coordinates
    Json::Value coords3d = root["atoms"]["coords"]["3d"];
    for(Json::ArrayIndex i = 0; i < coords3d.size(); i += 3){
        molecule->atom(i / 3)->setPosition(coords3d.get(i + 0, 0).asDouble(),
                                           coords3d.get(i + 1, 0).asDouble(),
                                           coords3d.get(i + 2, 0).asDouble());
    }

    // read bond connections
    Json::Value bondConnections = root["bonds"]["connections"];
    for(Json::ArrayIndex i = 0; i < bondConnections.size(); i += 2){
        molecule->addBond(bondConnections.get(i + 0, 0).asInt(),
                          bondConnections.get(i + 1, 0).asInt());
    }

    // read bond orders
    Json::Value bondOrders = root["bonds"]["order"];
    for(Json::ArrayIndex i = 0; i < bondOrders.size(); i++){
        molecule->bond(i)->setOrder(bondOrders.get(i, 1).asInt());
    }

    // read properties
    Json::Value properties = root["properties"];
    foreach(std::string propertyName, properties.getMemberNames()){
        Json::Value propertyValue = properties[propertyName];

        switch(propertyValue.type()){
        case Json::booleanValue:
            molecule->setData(propertyName, propertyValue.asBool());
            break;
        case Json::intValue:
            molecule->setData(propertyName, propertyValue.asInt());
            break;
        case Json::uintValue:
            molecule->setData(propertyName, propertyValue.asUInt());
            break;
        case Json::realValue:
            molecule->setData(propertyName, propertyValue.asDouble());
            break;
        case Json::stringValue:
            molecule->setData(propertyName, propertyValue.asString());
            break;
        default:
            break;
        }
    }

    // add molecule to file
    file->addMolecule(molecule);

    return true;
}

bool ChemJsonFileFormat::write(const chemkit::MoleculeFile *file, std::ostream &output)
{
    const boost::shared_ptr<chemkit::Molecule> &molecule = file->molecule();
    if(!molecule){
        setErrorString("File is empty");
        return false;
    }

    // start molecule block
    output << "{" << std::endl;

    // write version
    const int version = 0;
    output << "  \"chemical json\": " << version << "," << std::endl;

    // write molecule name
    std::string name = sanitizeJsonString(molecule->name());
    if(!name.empty()){
        output << "  \"name\": \"" << name << "\"," << std::endl;
    }

    // write molecular formula
    std::string formula = molecule->formula("spaced-formula");
    if(!formula.empty()){
        output << "  \"formula\": \"" << formula << "\"," << std::endl;
    }

    // write inchi formula
    std::string inchi = molecule->formula("inchi");
    if(!inchi.empty()){
        output << "  \"inchi\": \"" << inchi << "\"," << std::endl;
    }

    // start atom block
    output << "  \"atoms\": {" << std::endl;

    // write elements
    output << "    \"elements\": [";
    foreach(const chemkit::Atom *atom, molecule->atoms()){
        output << static_cast<int>(atom->atomicNumber());

        if(atom->index() != molecule->size() - 1){
            output << ", ";
        }
    }
    output << "]," << std::endl;

    // write coordinates
    output << "    \"coords\": {" << std::endl;
    output << "      \"3d\": [" << std::endl;
    foreach(const chemkit::Atom *atom, molecule->atoms()){
        output << "        " << atom->x() << ", " <<
                                atom->y() << ", " <<
                                atom->z();

        if(atom->index() != molecule->size() - 1){
            output << ",";
        }

        output << std::endl;
    }
    output << "      ]" << std::endl;
    output << "    }" << std::endl;

    // end atom block
    output << "  }," << std::endl;

    // start bond block
    output << "  \"bonds\": {" << std::endl;

    // write connections
    output << "    \"connections\": [" << std::endl;
    foreach(const chemkit::Bond *bond, molecule->bonds()){
        output << "      " << bond->atom1()->index() << ", " <<
                              bond->atom2()->index();

        if(bond->index() != molecule->bondCount() - 1){
            output << ", ";
        }

        output << std::endl;
    }
    output << "    ]," << std::endl;

    // write bond orders
    output << "    \"order\": [";
    foreach(const chemkit::Bond *bond, molecule->bonds()){
        output << static_cast<int>(bond->order());

        if(bond->index() != molecule->bondCount() - 1){
            output << ", ";
        }
    }
    output << "]" << std::endl;

    // end bond block
    output << "  }" << std::endl;

    // end molecule block
    output << "}" << std::endl;

    return true;
}
