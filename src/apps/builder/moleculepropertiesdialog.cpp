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

#include "moleculepropertiesdialog.h"
#include "ui_moleculepropertiesdialog.h"

#include <chemkit/molecule.h>

// --- Construction and Destruction ---------------------------------------- //
MoleculePropertiesDialog::MoleculePropertiesDialog(const chemkit::Molecule *molecule, QWidget *parent)
    : QDialog(parent),
      ui(new Ui::MoleculePropertiesDialog),
      m_molecule(molecule)
{
    ui->setupUi(this);

    if(!molecule)
        return;

    connect(ui->lineFormatComboBox, SIGNAL(currentIndexChanged(int)), SLOT(lineFormatChanged(int)));

    ui->nameValue->setText(molecule->name().c_str());
    ui->formulaValue->setText(formattedFormula(molecule));
    ui->atomCountValue->setText(QString::number(molecule->atomCount()));
    ui->bondCountValue->setText(QString::number(molecule->bondCount()));
    ui->molarMassValue->setText(QString("%1 g/mol").arg(molecule->mass()));
    ui->lineFormatValue->setText(molecule->formula("inchi").c_str());
}

MoleculePropertiesDialog::~MoleculePropertiesDialog()
{
    delete ui;
}

// --- Slots --------------------------------------------------------------- //
void MoleculePropertiesDialog::lineFormatChanged(int index)
{
    // inchi
    if(index == 0){
        ui->lineFormatValue->setText(m_molecule->formula("inchi").c_str());
    }

    // inchikey
    else if(index == 1){
        ui->lineFormatValue->setText(m_molecule->formula("inchikey").c_str());
    }

    // smiles
    else if(index == 2){
        ui->lineFormatValue->setText(m_molecule->formula("smiles").c_str());
    }
}

// --- Internal Methods ---------------------------------------------------- //
// Creates a HTML formula string for the molecule. The formatted
// string uses subscripts for the atom quantities.
QString MoleculePropertiesDialog::formattedFormula(const chemkit::Molecule *molecule) const
{
    QString formula;

    bool inSymbol = false;
    bool inNumber = false;

    foreach(const QChar &c, molecule->formula()){
        if(c.isLetter()){
            if(inSymbol){
                formula += c;
            }
            else{
                if(inNumber){
                    formula += "</sub>";
                    inNumber = false;
                }

                formula += c;
                inSymbol = true;
            }
        }
        else if(c.isNumber()){
            if(inNumber){
                formula += c;
            }
            else{
                if(inSymbol){
                    inSymbol = false;
                }

                formula += "<sub>";
                formula += c;
                inNumber = true;
            }
        }
        else{
            formula += c;
        }
    }

    return formula;
}
