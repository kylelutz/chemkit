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

#ifndef DISPLAYSETTINGSDOCK_H
#define DISPLAYSETTINGSDOCK_H

#include <QDockWidget>

class BuilderWindow;

namespace Ui
{
    class DisplaySettingsDock;
}

namespace chemkit
{
    class Molecule;
    class MoleculeWatcher;
    class GraphicsMoleculeItem;
}

class DisplaySettingsDock : public QDockWidget
{
    Q_OBJECT

    public:
        DisplaySettingsDock(BuilderWindow *builder);
        ~DisplaySettingsDock();

        void setShowHydrogens(bool showHydrogens);

    private slots:
        void moleculeDisplayTypeChanged(int index);
        void showHydrogensCheckClicked(bool checked);
        void showBondOrderCheckClicked(bool checked);
        void moleculeChanged(chemkit::Molecule *molecule);

    private:
        Ui::DisplaySettingsDock *ui;
        BuilderWindow *m_builder;
        chemkit::MoleculeWatcher *m_watcher;
        chemkit::GraphicsMoleculeItem *m_moleculeItem;
        bool m_showHydrogens;
};

#endif // DISPLAYSETTINGSDOCK_H
