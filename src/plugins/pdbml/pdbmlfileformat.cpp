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

#include "pdbmlfileformat.h"

#include <QtXml>

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

bool PdbmlFileFormat::read(QIODevice *iodev, chemkit::PolymerFile *file)
{
    QDomDocument doc;
    bool ok = doc.setContent(iodev, true);
    if(!ok){
        setErrorString(QString("PDBML parsing failed."));
        return false;
    }

    QDomElement datablockElement = doc.documentElement();

    chemkit::Polymer *polymer = new chemkit::Polymer;
    chemkit::PolymerChain *chain = 0;
    QHash<QString, chemkit::PolymerChain *> chainNameToChain;

    QString name = datablockElement.attribute("datablockName");
    if(!name.isEmpty()){
        polymer->setName(name);
    }

    QDomElement element = datablockElement.firstChildElement();
    while(!element.isNull()){
        // atoms
        if(element.tagName() == "atom_siteCategory"){
            QDomElement atomElement = element.firstChildElement();

            // residue and chain data
            chemkit::AminoAcid *residue = 0;
            int currentSequenceNumber = -1;
            QString currentChainName;

            while(!atomElement.isNull()){
                QString id = atomElement.attribute("id");
                Q_UNUSED(id);

                // read atom data
                QString symbol;
                QString group;
                QString x, y, z;
                QString chainName;
                int sequenceNumber = 0;
                QString atomType;
                QString residueSymbol;

                QDomElement atomDataElement = atomElement.firstChildElement();
                while(!atomDataElement.isNull()){
                    if(atomDataElement.tagName() == "type_symbol"){
                        symbol = atomDataElement.text();
                    }
                    else if(atomDataElement.tagName() == "Cartn_x"){
                        x = atomDataElement.text();
                    }
                    else if(atomDataElement.tagName() == "Cartn_y"){
                        y = atomDataElement.text();
                    }
                    else if(atomDataElement.tagName() == "Cartn_z"){
                        z = atomDataElement.text();
                    }
                    else if(atomDataElement.tagName() == "label_asym_id"){
                        chainName = atomDataElement.text();
                    }
                    else if(atomDataElement.tagName() == "label_seq_id"){
                        sequenceNumber = atomDataElement.text().toInt();
                    }
                    else if(atomDataElement.tagName() == "label_atom_id"){
                        atomType = atomDataElement.text();
                    }
                    else if(atomDataElement.tagName() == "label_comp_id"){
                        residueSymbol = atomDataElement.text();
                    }
                    else if(atomDataElement.tagName() == "group_PDB"){
                        group = atomDataElement.text();
                    }

                    atomDataElement = atomDataElement.nextSiblingElement();
                }

                // add atom and set its data
                chemkit::Atom *atom = polymer->addAtom(symbol);
                if(atom){
                    // atomic coordinates
                    atom->setPosition(x.toDouble(), y.toDouble(), z.toDouble());

                    if(group == "ATOM"){
                        // set residue
                        if(chainName != currentChainName){
                            chain = polymer->addChain();
                            currentChainName = chainName;
                            chainNameToChain[chainName] = chain;
                        }

                        if(sequenceNumber != currentSequenceNumber){
                            residue = new chemkit::AminoAcid(polymer);
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

                // go to next atom
                atomElement = atomElement.nextSiblingElement();
            }
        }

        // secondary structure
        else if(element.tagName() == "struct_confCategory"){
            QDomElement structElement = element.firstChildElement();

            while(!structElement.isNull()){
                QString chainName;
                int firstResidue = 0;
                int lastResidue = 0;
                chemkit::AminoAcid::Conformation conformation = chemkit::AminoAcid::Coil;

                QDomElement structDataElement = structElement.firstChildElement();
                while(!structDataElement.isNull()){
                    if(structDataElement.tagName() == "beg_label_seq_id"){
                        firstResidue = structDataElement.text().toInt();
                    }
                    else if(structDataElement.tagName() == "end_label_seq_id"){
                        lastResidue = structDataElement.text().toInt();
                    }
                    else if(structDataElement.tagName() == "beg_label_asym_id"){
                        chainName = structDataElement.text();
                    }
                    else if(structDataElement.tagName() == "conf_type_id"){
                        QString type = structDataElement.text();

                        if(type == "HELX_P"){
                            conformation = chemkit::AminoAcid::AlphaHelix;
                        }
                        else if(type == "TURN_P"){
                            conformation = chemkit::AminoAcid::BetaSheet;
                        }
                    }

                    structDataElement = structDataElement.nextSiblingElement();
                }

                chemkit::PolymerChain *chain = chainNameToChain.value(chainName, 0);
                if(chain){
                    for(int i = firstResidue; i < lastResidue; i++){
                        chemkit::AminoAcid *aminoAcid = static_cast<chemkit::AminoAcid *>(chain->residue(i - 1));
                        aminoAcid->setConformation(conformation);
                    }
                }

                structElement = structElement.nextSiblingElement();
            }
        }

        element = element.nextSiblingElement();
    }

    file->addPolymer(polymer);

    return true;
}
