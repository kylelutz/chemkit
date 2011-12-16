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

#include <QHash>
#include <QString>
#include <QVector>

#include <boost/algorithm/string.hpp>

#include <chemkit/atom.h>
#include <chemkit/polymer.h>
#include <chemkit/residue.h>
#include <chemkit/molecule.h>
#include <chemkit/aminoacid.h>
#include <chemkit/conformer.h>
#include <chemkit/nucleotide.h>
#include <chemkit/polymerfile.h>
#include <chemkit/polymerchain.h>

namespace {

// === PdbAtom ============================================================= //
class PdbAtom
{
public:
    PdbAtom(const char *data);

    int id;
    QString name;
    chemkit::Point3 position;
    int atomicNumber;
};

PdbAtom::PdbAtom(const char *data)
{
    // atom id
    sscanf(&data[5], "%d", &id);

    // atom name
    name.clear();
    for(int i = 13; i < 16; i++){
        name += data[i];
    }
    name = name.trimmed();

    // coordinates
    double x, y, z;
    sscanf(&data[31], "%lf%lf%lf", &x, &y, &z);
    position = chemkit::Point3(x, y, z);

    // atomic number
    QString symbol;
    for(int i = 77; i < 79 && isalpha(data[i]); i++){
        symbol += data[i];
    }
    atomicNumber = chemkit::Element::atomicNumber(symbol.trimmed().toStdString());
}

// === PdbResidue ========================================================== //
class PdbResidue
{
public:
    PdbResidue(const QString &name, int index);
    ~PdbResidue();

    void addAtom(PdbAtom *atom);
    QList<PdbAtom *> atoms() const;

    QString name() const;
    int index() const;

private:
    QString m_name;
    int m_index;
    QList<PdbAtom *> m_atoms;
};

PdbResidue::PdbResidue(const QString &name, int index)
    : m_name(name),
      m_index(index)
{
}

PdbResidue::~PdbResidue()
{
    qDeleteAll(m_atoms);
}

void PdbResidue::addAtom(PdbAtom *atom)
{
    m_atoms.append(atom);
}

QList<PdbAtom *> PdbResidue::atoms() const
{
    return m_atoms;
}

QString PdbResidue::name() const
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
    QString name() const;

    void addResidue(PdbResidue *residue);
    QList<PdbResidue *> residues() const;

    Type guessType() const;

private:
    char m_id;
    QString m_name;
    QList<PdbResidue *> m_residues;
};

PdbChain::PdbChain(char id)
    : m_id(id)
{
}

PdbChain::~PdbChain()
{
    qDeleteAll(m_residues);
}

char PdbChain::id() const
{
    return m_id;
}

QString PdbChain::name() const
{
    return m_name;
}

void PdbChain::addResidue(PdbResidue *residue)
{
    m_residues.append(residue);
}

QList<PdbResidue *> PdbChain::residues() const
{
    return m_residues;
}

PdbChain::Type PdbChain::guessType() const
{
    if(m_residues.isEmpty())
        return Protein;

    PdbResidue *residue = m_residues.first();

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
    QVector<chemkit::Point3> m_positions;
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

            m_positions.append(chemkit::Point3(x, y, z));
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

// === PdbFile ============================================================= //
class PdbFile
{
public:
    PdbFile();
    ~PdbFile();

    bool read(std::istream &input);

    void addChain(PdbChain *chain);
    void writePolymerFile(chemkit::PolymerFile *file);

private:
    QList<PdbChain *> m_chains;
    QList<PdbConformer *> m_conformers;
    QList<PdbConformation *> m_conformations;
};

PdbFile::PdbFile()
{
}


PdbFile::~PdbFile()
{
    qDeleteAll(m_chains);
    qDeleteAll(m_conformers);
    qDeleteAll(m_conformations);
}

bool PdbFile::read(std::istream &input)
{
    PdbChain *currentChain = 0;
    PdbResidue *currentResidue = 0;

    for(;;){
        char line[80];
        input.read(line, 79);
        if(input.eof()){
            break;
        }

        input.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

        if(strncmp("ATOM", line, 4) == 0){
            PdbAtom *atom = new PdbAtom(line);

            char chainId = line[21];
            if(!currentChain || currentChain->id() != chainId){
                currentChain = new PdbChain(chainId);
                addChain(currentChain);
            }

            int residueIndex = strtol(&line[22], 0, 10);
            if(!currentResidue || currentResidue->index() != residueIndex){
                QString name;
                for(int i = 17; i < 21; i++){
                    name += line[i];
                }

                name = name.trimmed();

                currentResidue = new PdbResidue(name, residueIndex);
                currentChain->addResidue(currentResidue);
            }

            currentResidue->addAtom(atom);
        }
        else if(strncmp("HELIX", line, 5) == 0 ||
                strncmp("SHEET", line, 5) == 0){
            m_conformations.append(new PdbConformation(line));
        }
        else if(strncmp("MODEL", line, 5) == 0 && !m_chains.isEmpty()){
            PdbConformer *conformer = new PdbConformer(input);
            m_conformers.append(conformer);
        }
    }

    return true;
}

void PdbFile::addChain(PdbChain *chain)
{
    m_chains.append(chain);
}

void PdbFile::writePolymerFile(chemkit::PolymerFile *file)
{
    if(m_chains.isEmpty()){
        return;
    }

    chemkit::Polymer *polymer = new chemkit::Polymer;

    QHash<int, chemkit::Atom *> atomIds;
    PdbChain::Type chainType = PdbChain::Protein;

    foreach(PdbChain *pdbChain, m_chains){
        chemkit::PolymerChain *chain = polymer->addChain();
        chainType = pdbChain->guessType();

        foreach(PdbResidue *pdbResidue, pdbChain->residues()){
            chemkit::AminoAcid *aminoAcid = 0;
            chemkit::Nucleotide *nucleotide = 0;
            chemkit::Residue *residue = 0;

            if(chainType == PdbChain::Protein){
                aminoAcid = new chemkit::AminoAcid(polymer);
                residue = aminoAcid;

                aminoAcid->setType(pdbResidue->name().toStdString());
            }
            else{
                nucleotide = new chemkit::Nucleotide(polymer);
                residue = nucleotide;

                QChar symbol;
                if(pdbResidue->name().length() == 1){
                    symbol = pdbResidue->name().at(0);
                    nucleotide->setSugarType(chemkit::Nucleotide::Ribose);
                }
                else if(pdbResidue->name().length() == 2 && pdbResidue->name().startsWith("D")){
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
                chemkit::Atom *atom = polymer->addAtom(pdbAtom->atomicNumber);
                if(!atom){
                    continue;
                }

                atomIds[pdbAtom->id] = atom;

                atom->setPosition(pdbAtom->position);
                residue->addAtom(atom);
                residue->setAtomType(atom, pdbAtom->name.toStdString());

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
        for(int i = 0; i < m_chains.size(); i++){
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
        chemkit::Conformer *conformer = polymer->addConformer();

        foreach(int atomId, atomIds.keys()){
            conformer->setPosition(atomIds[atomId], pdbConformer->position(atomId - 1));
        }
    }

    file->addPolymer(polymer);
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
