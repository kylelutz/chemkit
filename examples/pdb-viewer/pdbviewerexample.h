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

#ifndef PDBVIEWEREXAMPLE_H
#define PDBVIEWEREXAMPLE_H

#include <QtGui>

#include <chemkit/polymerfile.h>
#include <chemkit/graphicsview.h>
#include <chemkit/graphicsproteinitem.h>
#include <chemkit/graphicsmoleculeitem.h>
#include <chemkit/graphicsnucleicaciditem.h>

namespace Ui
{
    class PdbViewerWindow;
}

class PdbViewerWindow : public QMainWindow
{
    Q_OBJECT

    public:
        PdbViewerWindow(QWidget *parent = 0);
        ~PdbViewerWindow();

        void setFile(chemkit::PolymerFile *file);

    public slots:
        void openFile();
        void openFile(const QString &fileName);
        void quit();

    private:
        Ui::PdbViewerWindow *ui;
        chemkit::GraphicsView *m_view;
        chemkit::PolymerFile *m_file;
        chemkit::GraphicsProteinItem *m_proteinItem;
        chemkit::GraphicsNucleicAcidItem *m_nucleicAcidItem;
};

#endif // PDBVIEWEREXAMPLE_H
