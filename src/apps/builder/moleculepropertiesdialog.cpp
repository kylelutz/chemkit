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
    ui->lineFormatValue->setText(molecule->formula("inchi"));
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
        ui->lineFormatValue->setText(m_molecule->formula("inchi"));
    }

    // inchikey
    else if(index == 1){
        ui->lineFormatValue->setText(m_molecule->formula("inchikey"));
    }

    // smiles
    else if(index == 2){
        ui->lineFormatValue->setText(m_molecule->formula("smiles"));
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
