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
    void setMolecule(const boost::shared_ptr<chemkit::Molecule> &molecule);
    boost::shared_ptr<chemkit::Molecule> molecule() const;
    chemkit::GraphicsMoleculeItem* moleculeItem() const;
    chemkit::MoleculeEditor* editor() const;
    EnergyMinimizer* energyMinimizer() const;
    void beginMoleculeEdit();
    void endMoleculeEdit();

    // view
    chemkit::GraphicsView* view() const;

    // tools
    void setTool(const boost::shared_ptr<BuilderTool> &tool);
    boost::shared_ptr<BuilderTool> tool() const { return m_tool; }
    boost::shared_ptr<BuilderTool> navigateTool() const { return m_navigateTool; }
    boost::shared_ptr<BuilderTool> buildTool() const { return m_buildTool; }
    boost::shared_ptr<BuilderTool> manipulateTool() const { return m_manipulateTool; }
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
    boost::shared_ptr<chemkit::Molecule> m_molecule;
    chemkit::GraphicsMoleculeItem *m_moleculeItem;
    chemkit::MoleculeEditor *m_editor;
    boost::shared_ptr<BuilderTool> m_tool;
    boost::shared_ptr<BuilderTool> m_navigateTool;
    boost::shared_ptr<BuilderTool> m_buildTool;
    boost::shared_ptr<BuilderTool> m_manipulateTool;
    EnergyMinimizer *m_energyMinimizer;
    bool m_inMoleculeEdit;
};

#endif // BUILDERWINDOW_H
