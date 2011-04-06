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

#ifndef BUILDERWINDOW_H
#define BUILDERWINDOW_H

#include <QtGui>

#include <chemkit/chemkit.h>
#include <chemkit/molecule.h>
#include <chemkit/moleculeeditor.h>
#include <chemkit/moleculefileformat.h>
#include <chemkit/graphicsview.h>
#include <chemkit/graphicsmoleculeitem.h>

class BuilderTool;
class EnergyMinimizer;
class MoleculeLabelItem;

namespace Ui {
    class BuilderWindow;
}

class BuilderWindow : public QMainWindow
{
    Q_OBJECT

    public:
        // construction and destruction
        BuilderWindow(QWidget *parent = 0);
        ~BuilderWindow();

        // properties
        chemkit::MoleculeFile* file() const { return m_file; }

        // molecule
        void setMolecule(chemkit::Molecule *molecule);
        chemkit::Molecule* molecule() const;
        chemkit::GraphicsMoleculeItem* moleculeItem() const;
        chemkit::MoleculeEditor* editor() const;
        EnergyMinimizer* energyMinimizer() const;
        void beginMoleculeEdit();
        void endMoleculeEdit();

        // view
        chemkit::GraphicsView* view() const;

        // tools
        void setTool(BuilderTool *tool);
        BuilderTool* tool() const { return m_tool; }
        BuilderTool* navigateTool() const { return m_navigateTool; }
        BuilderTool* buildTool() const { return m_buildTool; }
        BuilderTool* manipulateTool() const { return m_manipulateTool; }
        void setCanCut(bool canCut);
        void setCanCopy(bool canCopy);
        void setCanPaste(bool canPaste);
        void setCanDelete(bool canDelete);

    public slots:
        void openFile(const QString &fileName);
        void openFile();
        void saveFile();
        void saveFileAs(const QString &fileName);
        void saveFileAs();
        void closeFile();
        void quit();
        void about();
        void undo();
        void redo();
        void cut();
        void copy();
        void paste();
        void del();
        void centerCamera();
        void setBackgroundColor(QAction *action);
        void setTool(QAction *action);
        void predictBonds();
        void adjustHydrogens();
        void moleculeProperties();

    signals:
        void fileChanged(chemkit::MoleculeFile *file);
        void moleculeChanged(chemkit::Molecule *molecule);
        void toolChanged(BuilderTool *tool);

    private slots:
        void minimizerStateChanged(int state);

    private:
        Ui::BuilderWindow *ui;
        chemkit::MoleculeFile *m_file;
        chemkit::Molecule *m_molecule;
        chemkit::GraphicsMoleculeItem *m_moleculeItem;
        chemkit::MoleculeEditor *m_editor;
        BuilderTool *m_tool;
        BuilderTool *m_navigateTool;
        BuilderTool *m_buildTool;
        BuilderTool *m_manipulateTool;
        EnergyMinimizer *m_energyMinimizer;
        bool m_inMoleculeEdit;
};

#endif // BUILDERWINDOW_H
