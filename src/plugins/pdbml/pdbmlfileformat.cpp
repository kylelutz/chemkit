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
        setErrorString("PDBML parsing failed.");
        return false;
    }

    QDomElement datablockElement = doc.documentElement();

    chemkit::Polymer *polymer = new chemkit::Polymer;
    chemkit::PolymerChain *chain = 0;
    QHash<QString, chemkit::PolymerChain *> chainNameToChain;

    QString name = datablockElement.attribute("datablockName");
    if(!name.isEmpty()){
        polymer->setName(name.toStdString());
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
                chemkit::Atom *atom = polymer->addAtom(symbol.toStdString());
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
                            residue->setType(residueSymbol.toStdString());
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
                        residue->setAtomType(atom, atomType.toStdString());
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
