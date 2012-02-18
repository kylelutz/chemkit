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

#include "pdbfileformat.h"

#include <boost/algorithm/string.hpp>

#include <chemkit/atom.h>
#include <chemkit/foreach.h>
#include <chemkit/polymer.h>
#include <chemkit/residue.h>
#include <chemkit/molecule.h>
#include <chemkit/aminoacid.h>
#include <chemkit/nucleotide.h>
#include <chemkit/polymerfile.h>
#include <chemkit/polymerchain.h>
#include <chemkit/coordinateset.h>
#include <chemkit/cartesiancoordinates.h>

namespace {

// === PdbAtom ============================================================= //
class PdbAtom
{
public:
    PdbAtom(const char *data);

    int id;
    std::string name;
    chemkit::Point3 position;
    chemkit::Element element;
};

PdbAtom::PdbAtom(const char *data)
{
    // atom id
    sscanf(&data[7], "%d", &id);

    // atom name
    name.clear();
    for(int i = 13; i < 16; i++){
        name += data[i];
    }
    boost::trim(name);

    // coordinates
    double x, y, z;
    sscanf(&data[31], "%lf%lf%lf", &x, &y, &z);
    position = chemkit::Point3(x, y, z);

    // atomic number
    std::string symbol;
    for(int i = 77; i < 79 && isalpha(data[i]); i++){
        symbol += tolower(data[i]);
    }
    boost::trim(symbol);
    symbol[0] = toupper(symbol[0]);
    element = chemkit::Element::fromSymbol(symbol);

    if(!element.isValid()){
        // try atomic number from name
        symbol = boost::to_lower_copy(name);
        symbol[0] = toupper(symbol[0]);
        element = chemkit::Element::fromSymbol(symbol);
    }
}

// === PdbResidue ========================================================== //
class PdbResidue
{
public:
    PdbResidue(const std::string &name, int index);
    ~PdbResidue();

    void addAtom(PdbAtom *atom);
    std::vector<PdbAtom *> atoms() const;

    std::string name() const;
    int index() const;

private:
    std::string m_name;
    int m_index;
    std::vector<PdbAtom *> m_atoms;
};

PdbResidue::PdbResidue(const std::string &name, int index)
    : m_name(name),
      m_index(index)
{
}

PdbResidue::~PdbResidue()
{
    foreach(PdbAtom *atom, m_atoms){
        delete atom;
    }
}

void PdbResidue::addAtom(PdbAtom *atom)
{
    m_atoms.push_back(atom);
}

std::vector<PdbAtom *> PdbResidue::atoms() const
{
    return m_atoms;
}

std::string PdbResidue::name() const
{
    return m_name;
}

int PdbResidue::index() const
{
    return m_index;
}

// === PdbChain ============================================================ //
class PdbChain
{
public:
    enum Type {
        Protein,
        DNA,
        RNA
    };

    PdbChain(char id);
    ~PdbChain();

    char id() const;
    std::string name() const;

    void addResidue(PdbResidue *residue);
    std::vector<PdbResidue *> residues() const;

    Type guessType() const;

private:
    char m_id;
    std::string m_name;
    std::vector<PdbResidue *> m_residues;
};

PdbChain::PdbChain(char id)
    : m_id(id)
{
}

PdbChain::~PdbChain()
{
    foreach(PdbResidue *residue, m_residues){
        delete residue;
    }
}

char PdbChain::id() const
{
    return m_id;
}

std::string PdbChain::name() const
{
    return m_name;
}

void PdbChain::addResidue(PdbResidue *residue)
{
    m_residues.push_back(residue);
}

std::vector<PdbResidue *> PdbChain::residues() const
{
    return m_residues;
}

PdbChain::Type PdbChain::guessType() const
{
    if(m_residues.empty())
        return Protein;

    PdbResidue *residue = m_residues.front();

    if(residue->name() == "DG" ||
       residue->name() == "DA" ||
       residue->name() == "DC" ||
       residue->name() == "DT"){
        return DNA;
    }

    return Protein;
}

// === PdbConformation ===================================================== //
class PdbConformation
{
public:
    PdbConformation(const char *data);

    chemkit::AminoAcid::Conformation type() const { return m_type; }
    char chain() const { return m_chain; }
    int firstResidue() const { return m_firstResidue; }
    int lastResidue() const { return m_lastResidue; }

private:
    chemkit::AminoAcid::Conformation m_type;
    char m_chain;
    int m_firstResidue;
    int m_lastResidue;
};

PdbConformation::PdbConformation(const char *data)
{
    m_type = chemkit::AminoAcid::Coil;

    if(strncmp("HELIX", data, 5) == 0){
        m_type = chemkit::AminoAcid::AlphaHelix;
        m_chain = data[19];
        sscanf(&data[21], "%d", &m_firstResidue);
        sscanf(&data[33], "%d", &m_lastResidue);
    }
    else if(strncmp("SHEET", data, 5) == 0){
        m_type = chemkit::AminoAcid::BetaSheet;
        m_chain = data[21];
        sscanf(&data[22], "%d", &m_firstResidue);
        sscanf(&data[33], "%d", &m_lastResidue);
    }
}

// === PdbConformer ======================================================== //
class PdbConformer
{
public:
    PdbConformer(std::istream &input);

    chemkit::Point3 position(int atom) const;

private:
    std::vector<chemkit::Point3> m_positions;
};

PdbConformer::PdbConformer(std::istream &input)
{
    for(;;){
        std::string line;
        std::getline(input, line);

        if(line.empty() || input.eof()){
            break;
        }

        if(boost::starts_with(line, "ATOM")){
            chemkit::Real x = boost::lexical_cast<chemkit::Real>(boost::trim_left_copy(line.substr(30, 8)));
            chemkit::Real y = boost::lexical_cast<chemkit::Real>(boost::trim_left_copy(line.substr(38, 8)));
            chemkit::Real z = boost::lexical_cast<chemkit::Real>(boost::trim_left_copy(line.substr(46, 8)));

            m_positions.push_back(chemkit::Point3(x, y, z));
        }
        else if(boost::starts_with(line, "ENDMDL")){
            break;
        }
    }
}

chemkit::Point3 PdbConformer::position(int atom) const
{
    return m_positions[atom];
}

// === PdbLigand =========================================================== //
class PdbLigand
{
public:
    PdbLigand(const std::string &name, int index);
    ~PdbLigand();

    std::string name() const;
    int index() const;
    void addAtom(PdbAtom *atom);
    std::vector<PdbAtom *> atoms() const;

private:
    int m_index;
    std::string m_name;
    std::vector<PdbAtom *> m_atoms;
};

PdbLigand::PdbLigand(const std::string &name, int index)
{
    m_name = name;
    m_index = index;
}

PdbLigand::~PdbLigand()
{
    foreach(PdbAtom *atom, m_atoms){
        delete atom;
    }
}

std::string PdbLigand::name() const
{
    return m_name;
}

int PdbLigand::index() const
{
    return m_index;
}

void PdbLigand::addAtom(PdbAtom *atom)
{
    m_atoms.push_back(atom);
}

std::vector<PdbAtom *> PdbLigand::atoms() const
{
    return m_atoms;
}

// === PdbFile ============================================================= //
class PdbFile
{
public:
    PdbFile();
    ~PdbFile();

    bool read(std::istream &input);

    void addChain(PdbChain *chain);
    void addLigand(PdbLigand *ligand);
    void addConnections(const std::vector<int> &connections);
    void writePolymerFile(chemkit::PolymerFile *file);

private:
    std::vector<PdbChain *> m_chains;
    std::vector<PdbConformer *> m_conformers;
    std::vector<PdbConformation *> m_conformations;
    std::vector<PdbLigand *> m_ligands;
    std::vector<std::vector<int> > m_connections;
    std::map<std::string, std::string> m_ligandNames;
    std::string m_title;
};

PdbFile::PdbFile()
{
}


PdbFile::~PdbFile()
{
    foreach(PdbChain *chain, m_chains)
        delete chain;
    foreach(PdbConformer *conformer, m_conformers)
        delete conformer;
    foreach(PdbConformation *conformation, m_conformations)
        delete conformation;
}

bool PdbFile::read(std::istream &input)
{
    PdbChain *currentChain = 0;
    PdbLigand *currentLigand = 0;
    PdbResidue *currentResidue = 0;

    for(;;){
        std::string lineString;
        std::getline(input, lineString);
        if(lineString.empty() || input.eof()){
            break;
        }

        const char *line = lineString.c_str();

        if(strncmp("ATOM", line, 4) == 0){
            PdbAtom *atom = new PdbAtom(line);

            char chainId = line[21];
            if(!currentChain || currentChain->id() != chainId){
                currentChain = new PdbChain(chainId);
                addChain(currentChain);
            }

            int residueIndex = strtol(&line[22], 0, 10);
            if(!currentResidue || currentResidue->index() != residueIndex){
                std::string name;
                for(int i = 17; i < 21; i++){
                    name += line[i];
                }

                boost::trim(name);

                currentResidue = new PdbResidue(name, residueIndex);
                currentChain->addResidue(currentResidue);
            }

            currentResidue->addAtom(atom);
        }
        else if(strncmp("HETATM", line, 6) == 0){
            PdbAtom *atom = new PdbAtom(line);

            int ligandIndex = strtol(&line[22], 0, 10);
            if(!currentLigand || currentLigand->index() != ligandIndex){
                std::string ligandName;
                for(int i = 17; i < 21; i++){
                    ligandName += line[i];
                }

                boost::trim(ligandName);

                currentLigand = new PdbLigand(ligandName, ligandIndex);
                addLigand(currentLigand);
            }

            currentLigand->addAtom(atom);
        }
        else if(strncmp("HELIX", line, 5) == 0 ||
                strncmp("SHEET", line, 5) == 0){
            m_conformations.push_back(new PdbConformation(line));
        }
        else if(strncmp("MODEL", line, 5) == 0 && !m_chains.empty()){
            PdbConformer *conformer = new PdbConformer(input);
            m_conformers.push_back(conformer);
        }
        else if(strncmp("CONECT", line, 6) == 0){
            std::string string = &line[7];
            boost::trim(string);

            std::vector<std::string> tokens;
            boost::split(tokens, string, boost::is_any_of(" "), boost::token_compress_on);

            std::vector<int> ids;

            foreach(const std::string &idString, tokens){
                try {
                    ids.push_back(boost::lexical_cast<int>(idString));
                }
                catch(boost::bad_lexical_cast&){
                }
            }

            addConnections(ids);
        }
        else if(strncmp("HETNAM", line, 6) == 0){
            std::string string = &line[7];
            boost::trim(string);

            std::vector<std::string> tokens;
            boost::split(tokens, string, boost::is_any_of(" "), boost::token_compress_on);

            if(tokens.size() > 1){
                std::string residueName = tokens[0];
                tokens.erase(tokens.begin());

                std::string name = boost::join(tokens, " ");

                m_ligandNames[residueName] = name;
            }
        }
        else if(strncmp("TITLE", line, 5) == 0){
            std::string title = &line[10];
            boost::trim_right(title);
            m_title += title;
        }
    }

    return true;
}

void PdbFile::addChain(PdbChain *chain)
{
    m_chains.push_back(chain);
}

void PdbFile::addLigand(PdbLigand *ligand)
{
    m_ligands.push_back(ligand);
}

void PdbFile::addConnections(const std::vector<int> &connections)
{
    m_connections.push_back(connections);
}

void PdbFile::writePolymerFile(chemkit::PolymerFile *file)
{
    boost::shared_ptr<chemkit::Polymer> polymer(new chemkit::Polymer);

    if(!m_title.empty()){
        polymer->setName(m_title);
    }

    std::map<int, chemkit::Atom *> atomIds;
    PdbChain::Type chainType = PdbChain::Protein;

    foreach(PdbChain *pdbChain, m_chains){
        chemkit::PolymerChain *chain = polymer->addChain();
        chainType = pdbChain->guessType();

        foreach(PdbResidue *pdbResidue, pdbChain->residues()){
            chemkit::AminoAcid *aminoAcid = 0;
            chemkit::Nucleotide *nucleotide = 0;
            chemkit::Residue *residue = 0;

            if(chainType == PdbChain::Protein){
                aminoAcid = new chemkit::AminoAcid(polymer.get());
                residue = aminoAcid;

                aminoAcid->setType(pdbResidue->name());
            }
            else{
                nucleotide = new chemkit::Nucleotide(polymer.get());
                residue = nucleotide;

                char symbol;
                if(pdbResidue->name().length() == 1){
                    symbol = pdbResidue->name().at(0);
                    nucleotide->setSugarType(chemkit::Nucleotide::Ribose);
                }
                else if(pdbResidue->name().length() == 2 && pdbResidue->name()[0] == 'D'){
                    symbol = pdbResidue->name().at(1);
                    nucleotide->setSugarType(chemkit::Nucleotide::Deoxyribose);
                }

                if(symbol == 'A'){
                    nucleotide->setType(chemkit::Nucleotide::Adenine);
                }
                else if(symbol == 'C'){
                    nucleotide->setType(chemkit::Nucleotide::Cytosine);
                }
                else if(symbol == 'G'){
                    nucleotide->setType(chemkit::Nucleotide::Guanine);
                }
                else if(symbol == 'T'){
                    nucleotide->setType(chemkit::Nucleotide::Thymine);
                }
                else if(symbol == 'U'){
                    nucleotide->setType(chemkit::Nucleotide::Uracil);
                }
            }

            foreach(PdbAtom *pdbAtom, pdbResidue->atoms()){
                chemkit::Atom *atom = polymer->addAtom(pdbAtom->element);
                if(!atom){
                    continue;
                }

                atomIds[pdbAtom->id] = atom;

                atom->setPosition(pdbAtom->position);
                residue->addAtom(atom);
                residue->setAtomType(atom, pdbAtom->name);

                if(chainType == PdbChain::Protein){
                    if(pdbAtom->name == "CA"){
                        aminoAcid->setAlphaCarbon(atom);
                    }
                    else if(pdbAtom->name == "N"){
                        aminoAcid->setAminoNitrogen(atom);
                    }
                    else if(pdbAtom->name == "C"){
                        aminoAcid->setCarbonylCarbon(atom);
                    }
                    else if(pdbAtom->name == "O"){
                        aminoAcid->setCarbonylOxygen(atom);
                    }
                }
            }

            chain->addResidue(residue);
        }
    }

    // set amino acid conformations (alpha helix or beta sheet)
    if(chainType == PdbChain::Protein){
        for(size_t i = 0; i < m_chains.size(); i++){
            PdbChain *pdbChain = m_chains[i];
            chemkit::PolymerChain *chain = polymer->chain(i);

            foreach(const PdbConformation *pdbConformation, m_conformations){
                if(pdbConformation->chain() == pdbChain->id()){
                    for(int residue = pdbConformation->firstResidue(); residue < pdbConformation->lastResidue(); residue++){
                        chemkit::AminoAcid *aminoAcid = static_cast<chemkit::AminoAcid *>(chain->residue(residue));

                        if(aminoAcid){
                            aminoAcid->setConformation(pdbConformation->type());
                        }
                    }
                }
            }
        }
    }

    // add conformers
    foreach(const PdbConformer *pdbConformer, m_conformers){
        chemkit::CartesianCoordinates *coordinates = new chemkit::CartesianCoordinates(polymer->size());

        for(size_t i = 0; i < polymer->size(); i++){
            coordinates->setPosition(i, pdbConformer->position(i));
        }

        polymer->addCoordinateSet(coordinates);
    }

    if(!polymer->isEmpty()){
        file->addPolymer(polymer);
    }

    // add ligands
    foreach(const PdbLigand *pdbLigand, m_ligands){
        boost::shared_ptr<chemkit::Molecule> ligand =
            boost::shared_ptr<chemkit::Molecule>(new chemkit::Molecule);

        std::map<std::string, std::string>::const_iterator iter = m_ligandNames.find(pdbLigand->name());
        if(iter != m_ligandNames.end()){
            ligand->setName(iter->second);
        }
        else{
            ligand->setName(pdbLigand->name());
        }

        foreach(const PdbAtom *pdbAtom, pdbLigand->atoms()){
            chemkit::Atom *atom = ligand->addAtom(pdbAtom->element);
            if(!atom){
                continue;
            }

            atom->setPosition(pdbAtom->position);

            atomIds[pdbAtom->id] = atom;
        }

        file->addLigand(ligand);
    }

    // add connections
    foreach(const std::vector<int> &connections, m_connections){
        if(connections.size() < 2){
            continue;
        }

        int atomIdA = connections[0];
        chemkit::Atom *a = atomIds[atomIdA];
        if(!a){
            continue;
        }

        chemkit::Molecule *molecule = a->molecule();

        for(size_t i = 1; i < connections.size(); i++){
            int atomIdB = connections[i];
            chemkit::Atom *b = atomIds[atomIdB];
            if(!b){
                continue;
            }

            molecule->addBond(a, b);
        }
    }
}

} // end anonymous namespace

PdbFileFormat::PdbFileFormat()
    : chemkit::PolymerFileFormat("pdb")
{
}

bool PdbFileFormat::read(std::istream &input, chemkit::PolymerFile *file)
{
    PdbFile pdb;
    bool ok = pdb.read(input);
    if(!ok){
        return false;
    }

    pdb.writePolymerFile(file);
    return true;
}
