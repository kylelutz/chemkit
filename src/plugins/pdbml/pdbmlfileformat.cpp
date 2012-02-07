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

#include "pdbmlfileformat.h"

#include "../../3rdparty/rapidxml/rapidxml.hpp"

#include <chemkit/atom.h>
#include <chemkit/polymer.h>
#include <chemkit/aminoacid.h>
#include <chemkit/polymerfile.h>
#include <chemkit/polymerchain.h>

PdbmlFileFormat::PdbmlFileFormat()
    : chemkit::PolymerFileFormat("pdbml")
{
}

PdbmlFileFormat::~PdbmlFileFormat()
{
}

bool PdbmlFileFormat::read(std::istream &input, chemkit::PolymerFile *file)
{
    // read file data into a string
    std::string data((std::istreambuf_iterator<char>(input)),
                     std::istreambuf_iterator<char>());

    // parse document
    rapidxml::xml_document<> doc;
    doc.parse<0>(const_cast<char *>(data.c_str()));

    // parse polymers
    boost::shared_ptr<chemkit::Polymer> polymer(new chemkit::Polymer);
    chemkit::PolymerChain *chain = 0;
    std::map<std::string, chemkit::PolymerChain *> nameToChain;

    rapidxml::xml_node<> *datablockNode = doc.first_node("PDBx:datablock");
    rapidxml::xml_attribute<> *nameAttribute = datablockNode->first_attribute("datablockName");
    if(nameAttribute && nameAttribute->value()){
        polymer->setName(nameAttribute->value());
    }

    rapidxml::xml_node<> *node = datablockNode->first_node();
    while(node){
        // atoms
        if(strcmp(node->name(), "PDBx:atom_siteCategory") == 0){
            rapidxml::xml_node<> *atomNode = node->first_node("PDBx:atom_site");

            // residue and chain data
            chemkit::AminoAcid *residue = 0;
            int currentSequenceNumber = -1;
            std::string currentChainName;

            while(atomNode){
                rapidxml::xml_attribute<> *atomIdAttribute = atomNode->first_attribute("id");

                // read id
                std::string id;
                if(atomIdAttribute && atomIdAttribute->value()){
                    id = atomIdAttribute->value();
                }
                CHEMKIT_UNUSED(id);

                // read atom data
                std::string symbol;
                std::string group;
                std::string x, y, z;
                std::string chainName;
                int sequenceNumber = 0;
                std::string atomType;
                std::string residueSymbol;

                rapidxml::xml_node<> *atomDataNode = atomNode->first_node();
                while(atomDataNode){
                    if(strcmp(atomDataNode->name(), "PDBx:type_symbol") == 0){
                        symbol = atomDataNode->value();
                    }
                    else if(strcmp(atomDataNode->name(), "PDBx:Cartn_x") == 0){
                        x = atomDataNode->value();
                    }
                    else if(strcmp(atomDataNode->name(), "PDBx:Cartn_y") == 0){
                        y = atomDataNode->value();
                    }
                    else if(strcmp(atomDataNode->name(), "PDBx:Cartn_z") == 0){
                        z = atomDataNode->value();
                    }
                    else if(strcmp(atomDataNode->name(), "PDBx:label_asym_id") == 0){
                        chainName = atomDataNode->value();
                    }
                    else if(strcmp(atomDataNode->name(), "PDBx:label_seq_id") == 0){
                        if(atomDataNode->value() && atomDataNode->value_size()){
                            sequenceNumber = boost::lexical_cast<int>(atomDataNode->value());
                        }
                    }
                    else if(strcmp(atomDataNode->name(), "PDBx:label_atom_id") == 0){
                        atomType = atomDataNode->value();
                    }
                    else if(strcmp(atomDataNode->name(), "PDBx:label_comp_id") == 0){
                        residueSymbol = atomDataNode->value();
                    }
                    else if(strcmp(atomDataNode->name(), "PDBx:group_PDB") == 0){
                        group = atomDataNode->value();
                    }

                    atomDataNode = atomDataNode->next_sibling();
                }

                // add atom and set its data
                chemkit::Atom *atom = polymer->addAtom(symbol);
                if(atom){
                    // atomic coordinates
                    if(!x.empty() && !y.empty() && !z.empty()){
                        atom->setPosition(boost::lexical_cast<chemkit::Real>(x),
                                          boost::lexical_cast<chemkit::Real>(y),
                                          boost::lexical_cast<chemkit::Real>(z));
                    }

                    if(group == "ATOM"){
                        // set residue
                        if(chainName != currentChainName){
                            chain = polymer->addChain();
                            currentChainName = chainName;
                            nameToChain[chainName] = chain;
                        }

                        if(sequenceNumber != currentSequenceNumber){
                            residue = new chemkit::AminoAcid(polymer.get());
                            residue->setType(residueSymbol);
                            chain->addResidue(residue);
                            currentSequenceNumber = sequenceNumber;
                        }

                        if(atomType == "CA"){
                            residue->setAlphaCarbon(atom);
                        }
                        else if(atomType == "C"){
                            residue->setCarbonylCarbon(atom);
                        }
                        else if(atomType == "O"){
                            residue->setCarbonylOxygen(atom);
                        }
                        else if(atomType == "N"){
                            residue->setAminoNitrogen(atom);
                        }
                        residue->setAtomType(atom, atomType);
                    }
                }

                // go to next atom node
                atomNode = atomNode->next_sibling("PDBx:atom_site");
            }
        }

        // secondary structure
        else if(strcmp(node->name(), "PDBx:struct_confCategory") == 0){
            rapidxml::xml_node<> *structNode = node->first_node();

            while(structNode){
                std::string chainName;
                int firstResidue = 0;
                int lastResidue = 0;
                chemkit::AminoAcid::Conformation conformation = chemkit::AminoAcid::Coil;

                rapidxml::xml_node<> *structDataNode = structNode->first_node();
                while(structDataNode){
                    if(strcmp(structDataNode->name(), "PDBx:beg_label_seq_id") == 0){
                        if(structDataNode->value()){
                            firstResidue = boost::lexical_cast<int>(structDataNode->value());
                        }
                    }
                    else if(strcmp(structDataNode->name(), "PDBx:end_label_seq_id") == 0){
                        if(structDataNode->value()){
                            lastResidue = boost::lexical_cast<int>(structDataNode->value());
                        }
                    }
                    else if(strcmp(structDataNode->name(), "PDBx:beg_label_asym_id") == 0){
                        chainName = structDataNode->value();
                    }
                    else if(strcmp(structDataNode->name(), "PDBx:conf_type_id") == 0){
                        std::string type;

                        if(structDataNode->value()){
                            type = structDataNode->value();
                        }

                        if(type == "HELX_P"){
                            conformation = chemkit::AminoAcid::AlphaHelix;
                        }
                        else if(type == "TURN_P"){
                            conformation = chemkit::AminoAcid::BetaSheet;
                        }
                    }

                    structDataNode = structDataNode->next_sibling();
                }

                std::map<std::string, chemkit::PolymerChain *>::const_iterator iter = nameToChain.find(chainName);
                if(iter != nameToChain.end()){
                    chemkit::PolymerChain *chain = iter->second;

                    for(int i = firstResidue; i < lastResidue; i++){
                        chemkit::AminoAcid *aminoAcid = static_cast<chemkit::AminoAcid *>(chain->residue(i - 1));
                        aminoAcid->setConformation(conformation);
                    }
                }

                structNode = structNode->next_sibling();
            }
        }

        node = node->next_sibling();
    }

    file->addPolymer(polymer);

    return true;
}
