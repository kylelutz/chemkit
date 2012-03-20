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

#include "cmlfileformat.h"

#include "../../3rdparty/rapidxml/rapidxml.hpp"

#include <boost/make_shared.hpp>

#include <chemkit/atom.h>
#include <chemkit/bond.h>
#include <chemkit/foreach.h>
#include <chemkit/molecule.h>
#include <chemkit/moleculefile.h>
#include <chemkit/diagramcoordinates.h>
#include <chemkit/cartesiancoordinates.h>

CmlFileFormat::CmlFileFormat()
    : chemkit::MoleculeFileFormat("cml")
{
}

CmlFileFormat::~CmlFileFormat()
{
}

bool CmlFileFormat::read(std::istream &input, chemkit::MoleculeFile *file)
{
    // read file data into a string
    std::string data((std::istreambuf_iterator<char>(input)),
                      std::istreambuf_iterator<char>());

    // parse document
    rapidxml::xml_document<> doc;
    try {
        doc.parse<0>(const_cast<char *>(data.c_str()));
    }
    catch(rapidxml::parse_error &e){
        setErrorString(std::string("XML parse error: ") + e.what());
        return false;
    }

    // parse molecules
    boost::shared_ptr<chemkit::Molecule> molecule;
    rapidxml::xml_node<> *moleculeNode = doc.first_node("molecule");
    while(moleculeNode){
        molecule = boost::make_shared<chemkit::Molecule>();

        chemkit::DiagramCoordinates *diagramCoordinates = 0;
        chemkit::CartesianCoordinates *cartesianCoordinates = 0;

        // parse name
        rapidxml::xml_node<> *nameNode = moleculeNode->first_node("name");
        if(nameNode && nameNode->value()){
            molecule->setName(nameNode->value());
        }

        // parse atoms
        rapidxml::xml_node<> *atomArrayNode = moleculeNode->first_node("atomArray");
        if(atomArrayNode){
            rapidxml::xml_node<> *atomNode = atomArrayNode->first_node("atom");

            while(atomNode){
                chemkit::Point2f point2(0, 0);
                chemkit::Point3 point3(0, 0, 0);

                rapidxml::xml_attribute<> *attr = atomNode->first_attribute();
                while(attr){
                    if(strcmp(attr->name(), "elementType") == 0){
                        molecule->addAtom(attr->value());
                    }
                    else if(strcmp(attr->name(), "x2") == 0){
                        point2[0] = static_cast<float>(strtod(attr->value(), 0));
                    }
                    else if(strcmp(attr->name(), "y2") == 0){
                        point2[1] = static_cast<float>(strtod(attr->value(), 0));
                    }
                    else if(strcmp(attr->name(), "x3") == 0){
                        point3[0] = strtod(attr->value(), 0);
                    }
                    else if(strcmp(attr->name(), "y3") == 0){
                        point3[1] = strtod(attr->value(), 0);
                    }
                    else if(strcmp(attr->name(), "z3") == 0){
                        point3[2] = strtod(attr->value(), 0);
                    }

                    attr = attr->next_attribute();
                }

                if(cartesianCoordinates){
                    cartesianCoordinates->append(point3);
                }
                else if(!point3.isZero()){
                    cartesianCoordinates = new chemkit::CartesianCoordinates(molecule->size() - 1);
                    cartesianCoordinates->append(point3);
                }

                if(diagramCoordinates){
                    diagramCoordinates->append(point2);
                }
                else if(!point2.isZero()){
                    diagramCoordinates = new chemkit::DiagramCoordinates(molecule->size() - 1);
                    diagramCoordinates->append(point2);
                }

                atomNode = atomNode->next_sibling("atom");
            }
        }

        // parse bonds
        rapidxml::xml_node<> *bondArrayNode = moleculeNode->first_node("bondArray");
        if(bondArrayNode){
            rapidxml::xml_node<> *bondNode = bondArrayNode->first_node("bond");
            while(bondNode){
                rapidxml::xml_attribute<> *atomRefs2Attr = bondNode->first_attribute("atomRefs2");
                if(atomRefs2Attr && atomRefs2Attr->value()){
                    unsigned int atom1;
                    unsigned int atom2;
                    int count = sscanf(atomRefs2Attr->value(), " %*c%u %*c%u", &atom1, &atom2);
                    if(count == 2){
                        rapidxml::xml_attribute<> *orderAttr = bondNode->first_attribute("order");
                        chemkit::Bond::BondOrderType bondOrder = chemkit::Bond::Single;

                        if(orderAttr && orderAttr->value()){
                            bondOrder = strtol(orderAttr->value(), 0, 10);
                        }

                        molecule->addBond(molecule->atom(atom1 - 1),
                                          molecule->atom(atom2 - 1),
                                          bondOrder);
                    }
                }

                bondNode = bondNode->next_sibling("bond");
            }
        }

        // add coordinate sets
        if(cartesianCoordinates){
            molecule->addCoordinateSet(cartesianCoordinates);
            cartesianCoordinates = 0;
        }

        if(diagramCoordinates){
            molecule->addCoordinateSet(diagramCoordinates);
            diagramCoordinates = 0;
        }

        // add molecule to file
        file->addMolecule(molecule);
        molecule.reset();

        // move to next molecule
        moleculeNode = moleculeNode->next_sibling("molecule");
    }

    return true;
}

bool CmlFileFormat::write(const chemkit::MoleculeFile *file, std::ostream &output)
{
    output << "<?xml version=\"1.0\"?>\n";

    // write each molecule
    foreach(const boost::shared_ptr<chemkit::Molecule> &molecule, file->molecules()){
        output << "<molecule>\n";

        // write molecule name
        if(!molecule->name().empty()){
            output << "  <name>" << molecule->name() << "</name>\n";
        }

        // write atom array
        if(molecule->atomCount() != 0){
            output << "  <atomArray>\n";

            foreach(const chemkit::Atom *atom, molecule->atoms()){
                output << "    <atom id=\"a" << atom->index() + 1 << "\""
                       << " elementType=\"" << atom->symbol() << "\""
                       << " x3=\"" << atom->x() << "\""
                       << " y3=\"" << atom->y() << "\""
                       << " z3=\"" << atom->z() << "\""
                       << "/>\n";
            }

            output << "  </atomArray>\n";
        }

        // write bond array
        if(molecule->bondCount() != 0){
            output << "  <bondArray>\n";

            foreach(const chemkit::Bond *bond, molecule->bonds()){
                output << "    <bond atomRefs2=\""
                       << "a" << bond->atom1()->index() + 1 << " "
                       << "a" << bond->atom2()->index() + 1 << "\""
                       << " order=\"" << static_cast<int>(bond->order())
                       << "\"/>\n";
            }

            output << "  </bondArray>\n";
        }

        output << "</molecule>\n";
    }

    return true;
}
