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

#ifndef MOLECULEPROPERTIESDIALOG_H
#define MOLECULEPROPERTIESDIALOG_H

#include <QDialog>

namespace Ui
{
    class MoleculePropertiesDialog;
}

namespace chemkit
{
    class Molecule;
}

class MoleculePropertiesDialog : public QDialog
{
    Q_OBJECT

    public:
        explicit MoleculePropertiesDialog(const chemkit::Molecule *molecule, QWidget *parent = 0);
        ~MoleculePropertiesDialog();

    private slots:
        void lineFormatChanged(int index);

    private:
        QString formattedFormula(const chemkit::Molecule *molecule) const;

    private:
        Ui::MoleculePropertiesDialog *ui;
        const chemkit::Molecule *m_molecule;
};

#endif // MOLECULEPROPERTIESDIALOG_H
