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

#include "builderwindow.h"
#include "ui_builderwindow.h"

#include <chemkit/bondpredictor.h>
#include <chemkit/graphicscamera.h>

#include "buildertool.h"
#include "navigatetool.h"
#include "buildtool.h"
#include "manipulatetool.h"
#include "energyminimizer.h"
#include "moleculelistdock.h"
#include "toolsettingsdock.h"
#include "displaysettingsdock.h"
#include "energyminimizationdock.h"
#include "moleculepropertiesdialog.h"

// --- Construction and Destruction ---------------------------------------- //
BuilderWindow::BuilderWindow(QWidget *parent)
    : QMainWindow(parent),
      ui(new Ui::BuilderWindow)
{
    // properties
    m_file = 0;
    m_molecule = 0;
    m_inMoleculeEdit = false;

    // setup ui
    ui->setupUi(this);
    connect(ui->actionOpen, SIGNAL(triggered()), SLOT(openFile()));
    connect(ui->actionSave, SIGNAL(triggered()), SLOT(saveFile()));
    connect(ui->actionSaveAs, SIGNAL(triggered()), SLOT(saveFileAs()));
    connect(ui->actionClose, SIGNAL(triggered()), SLOT(closeFile()));
    connect(ui->actionQuit, SIGNAL(triggered()), SLOT(quit()));
    connect(ui->actionAbout, SIGNAL(triggered()), SLOT(about()));
    connect(ui->actionUndo, SIGNAL(triggered()), SLOT(undo()));
    connect(ui->actionRedo, SIGNAL(triggered()), SLOT(redo()));
    connect(ui->actionCut, SIGNAL(triggered()), SLOT(cut()));
    connect(ui->actionCopy, SIGNAL(triggered()), SLOT(copy()));
    connect(ui->actionPaste, SIGNAL(triggered()), SLOT(paste()));
    connect(ui->actionDelete, SIGNAL(triggered()), SLOT(del()));
    connect(ui->actionCenterCamera, SIGNAL(triggered()), SLOT(centerCamera()));
    connect(ui->actionPredictBonds, SIGNAL(triggered()), SLOT(predictBonds()));
    connect(ui->actionAdjustHydrogens, SIGNAL(triggered()), SLOT(adjustHydrogens()));
    connect(ui->actionMoleculeProperties, SIGNAL(triggered()), SLOT(moleculeProperties()));

    QActionGroup *toolGroup = new QActionGroup(this);
    toolGroup->addAction(ui->actionNavigate);
    toolGroup->addAction(ui->actionBuild);
    toolGroup->addAction(ui->actionManipulate);
    toolGroup->setExclusive(true);
    connect(toolGroup, SIGNAL(triggered(QAction*)), SLOT(setTool(QAction*)));

    QActionGroup *backgroundColorGroup = new QActionGroup(this);
    backgroundColorGroup->addAction(ui->actionBackgroundBlack);
    backgroundColorGroup->addAction(ui->actionBackgroundWhite);
    backgroundColorGroup->addAction(ui->actionBackgroundGray);
    backgroundColorGroup->addAction(ui->actionBackgroundOther);
    backgroundColorGroup->setExclusive(true);
    connect(backgroundColorGroup, SIGNAL(triggered(QAction*)), SLOT(setBackgroundColor(QAction*)));

    // setup icons for menus and toolbars
    ui->actionOpen->setIcon(QIcon::fromTheme("document-open", style()->standardIcon(QStyle::SP_DialogOpenButton)));
    ui->actionSave->setIcon(QIcon::fromTheme("document-save", style()->standardIcon(QStyle::SP_DialogSaveButton)));
    ui->actionSaveAs->setIcon(QIcon::fromTheme("document-save-as"));
    ui->actionClose->setIcon(QIcon::fromTheme("window-close", style()->standardIcon(QStyle::SP_DirClosedIcon)));
    ui->actionQuit->setIcon(QIcon::fromTheme("application-exit", style()->standardIcon(QStyle::SP_DialogCloseButton)));
    ui->actionUndo->setIcon(QIcon::fromTheme("edit-undo", style()->standardIcon(QStyle::SP_ArrowBack)));
    ui->actionRedo->setIcon(QIcon::fromTheme("edit-redo", style()->standardIcon(QStyle::SP_ArrowForward)));
    ui->actionCut->setIcon(QIcon::fromTheme("edit-cut"));
    ui->actionCopy->setIcon(QIcon::fromTheme("edit-copy"));
    ui->actionPaste->setIcon(QIcon::fromTheme("edit-paste"));
    ui->actionDelete->setIcon(QIcon::fromTheme("edit-delete"));
    ui->actionAbout->setIcon(QIcon::fromTheme("help-about"));

    // setup view
    m_moleculeItem = 0;

    // setup molecule editor
    m_editor = new chemkit::MoleculeEditor;
    ui->actionUndo->setEnabled(false);
    ui->actionRedo->setEnabled(false);
    ui->actionCut->setEnabled(false);
    ui->actionCopy->setEnabled(false);
    ui->actionPaste->setEnabled(false);
    ui->actionDelete->setEnabled(false);
    connect(m_editor, SIGNAL(canUndoChanged(bool)), ui->actionUndo, SLOT(setEnabled(bool)));
    connect(m_editor, SIGNAL(canRedoChanged(bool)), ui->actionRedo, SLOT(setEnabled(bool)));
    connect(m_editor, SIGNAL(canPasteChanged(bool)), ui->actionPaste, SLOT(setEnabled(bool)));

    // setup energy minimizer
    m_energyMinimizer = new EnergyMinimizer;
    connect(m_energyMinimizer, SIGNAL(stateChanged(int)), SLOT(minimizerStateChanged(int)));

    // central view widget
    m_view = new chemkit::GraphicsView;
    setCentralWidget(m_view);

    // setup tools
    m_tool = 0;
    m_navigateTool = new NavigateTool(this);
    m_buildTool = new BuildTool(this);
    m_manipulateTool = new ManipulateTool(this);
    setTool(m_navigateTool);

    // dock widgets
    QDockWidget *dockWidget;

    // tool settings dock widget
    dockWidget = new ToolSettingsDock(this);
    addDockWidget(Qt::LeftDockWidgetArea, dockWidget);
    ui->menuView->addAction(dockWidget->toggleViewAction());

    // display settings dock widget
    dockWidget = new DisplaySettingsDock(this);
    addDockWidget(Qt::LeftDockWidgetArea, dockWidget);
    ui->menuView->addAction(dockWidget->toggleViewAction());

    // energy minimization dock widget
    dockWidget = new EnergyMinimizationDock(this);
    addDockWidget(Qt::LeftDockWidgetArea, dockWidget);
    ui->menuView->addAction(dockWidget->toggleViewAction());

    // molecule list dock widget
    dockWidget = new MoleculeListDock(this);
    addDockWidget(Qt::RightDockWidgetArea, dockWidget);
    dockWidget->setVisible(false);
    ui->menuView->addAction(dockWidget->toggleViewAction());

    //setMolecule(0);
    setMolecule(new chemkit::Molecule);

    setTool(m_buildTool);
}

BuilderWindow::~BuilderWindow()
{
    m_view->setTool(0);
    delete m_navigateTool;
    delete m_buildTool;
    delete m_manipulateTool;

//    delete m_editor;
//    delete m_molecule;
//    delete m_energyMinimizer;

    delete m_view;
    delete ui;
}

void BuilderWindow::setTool(BuilderTool *tool)
{
    if(tool == m_tool){
        return; // no change
    }

    m_tool = tool;
    m_view->setTool(tool);
    emit toolChanged(tool);

    if(tool == m_navigateTool)
        ui->actionNavigate->setChecked(true);
    else if(tool == m_buildTool)
        ui->actionBuild->setChecked(true);
    else if(tool == m_manipulateTool)
        ui->actionManipulate->setChecked(true);
}

void BuilderWindow::setMolecule(chemkit::Molecule *molecule)
{
    if(m_molecule == molecule){
        return;
    }

    // set molecule
    m_molecule = molecule;

    // remove old molecule item
    if(m_moleculeItem){
        m_view->deleteItem(m_moleculeItem);
        m_moleculeItem = 0;
    }

    // add new molecule item
    if(molecule){
        m_moleculeItem = new chemkit::GraphicsMoleculeItem(molecule);
        m_view->addItem(m_moleculeItem);
    }

    // reset editor
    m_editor->setMolecule(molecule);

    // reset energy minimizer
    m_energyMinimizer->setMolecule(molecule);

    centerCamera();

    // notify observers
    emit moleculeChanged(molecule);
}

chemkit::Molecule* BuilderWindow::molecule() const
{
    return m_molecule;
}

chemkit::GraphicsMoleculeItem* BuilderWindow::moleculeItem() const
{
    return m_moleculeItem;
}

chemkit::MoleculeEditor* BuilderWindow::editor() const
{
    return m_editor;
}

EnergyMinimizer* BuilderWindow::energyMinimizer() const
{
    return m_energyMinimizer;
}

void BuilderWindow::beginMoleculeEdit()
{
    m_editor->beginEdit();
    m_inMoleculeEdit = true;
}

void BuilderWindow::endMoleculeEdit()
{
    m_editor->endEdit();
    m_inMoleculeEdit = false;

    if(m_energyMinimizer->state() == EnergyMinimizer::Running ||
       m_energyMinimizer->state() == EnergyMinimizer::Converged ||
       m_energyMinimizer->state() == EnergyMinimizer::UpdateReady){
        m_energyMinimizer->start();
    }
}

// --- Slots --------------------------------------------------------------- //
void BuilderWindow::openFile(const QString &fileName)
{
    // close current file
    closeFile();

    // open and read file
    chemkit::ChemicalFile *file = new chemkit::ChemicalFile(fileName);
    if(!file->read()){
        QMessageBox::critical(this, "Error", QString("Error opening file: %1").arg(file->errorString()));
        delete file;
        return;
    }

    // set new file
    m_file = file;
    emit fileChanged(m_file);

    // set new molecule
    if(file->moleculeCount() > 0)
        setMolecule(file->molecule());

    if(!m_molecule || m_molecule->isEmpty()){
        setTool(m_buildTool);
    }
    else{
        setTool(m_navigateTool);
    }
}

void BuilderWindow::openFile()
{
    QStringList formatList = chemkit::ChemicalFile::formats();
    formatList.sort();

    QString formats;
    foreach(const QString &format, formatList)
        formats += QString("*.%1 ").arg(format);

    QString fileName = QFileDialog::getOpenFileName(this, tr("Open File"), 0, QString("Molecule Files (%1);;All Files (*.*)").arg(formats));

    if(!fileName.isEmpty()){
        openFile(fileName);
    }
}

void BuilderWindow::saveFile()
{
    if(!m_file || m_file->fileName().isEmpty()){
        saveFileAs();
        return;
    }

    bool ok = m_file->write();

    if(!ok){
        QMessageBox::critical(this, "Error", QString("Error saving file: %1").arg(m_file->errorString()));
    }
}

void BuilderWindow::saveFileAs(const QString &fileName)
{
    if(!m_file){
        m_file = new chemkit::ChemicalFile(fileName);
        m_file->addMolecule(m_molecule);
    }

    m_file->setFileName(fileName);
    saveFile();
}

void BuilderWindow::saveFileAs()
{
    QString fileName = QFileDialog::getSaveFileName(this, tr("Save File As"));

    if(!fileName.isEmpty())
        saveFileAs(fileName);
}

void BuilderWindow::closeFile()
{
    // remove molecule
    setMolecule(0);

    // remove file
    delete m_file;
    m_file = 0;
    emit fileChanged(m_file);
}

void BuilderWindow::quit()
{
    closeFile();
    qApp->quit();
}

void BuilderWindow::about()
{
    QString text;
    text += "<h2>chemkit-builder</h2>";
    text += "A molecular editor built ";
    text += "with the chemkit library. ";
    text += "See <a href=http://www.chemkit.org>http://www.chemkit.org</a> for more information.";

    QMessageBox::about(this, "About", text);
}

void BuilderWindow::undo()
{
    m_editor->undo();
    view()->update();
}

void BuilderWindow::redo()
{
    m_editor->redo();
    view()->update();
}

void BuilderWindow::cut()
{
    m_tool->cut();
}

void BuilderWindow::copy()
{
    m_tool->copy();
}

void BuilderWindow::paste()
{
    m_tool->paste();
}

void BuilderWindow::del()
{
    m_tool->del();
}

void BuilderWindow::minimizerStateChanged(int state)
{
    if(m_inMoleculeEdit){
        return;
    }

    if(state == EnergyMinimizer::UpdateReady){

        m_editor->beginEdit();
        // update atom positions
        chemkit::ForceField *forceField = m_energyMinimizer->forceField();
        foreach(chemkit::ForceFieldAtom *forceFieldAtom, forceField->atoms()){
            m_editor->setAtomPosition(const_cast<chemkit::Atom *>(forceFieldAtom->atom()), forceFieldAtom->position());
        }
        m_editor->endEdit();

        // run another step
        m_energyMinimizer->start();
    }
}

void BuilderWindow::centerCamera()
{
    if(m_molecule){
        view()->camera()->lookAt(m_molecule->center());
        view()->update();
    }
}

void BuilderWindow::setBackgroundColor(QAction *action)
{
    if(action == ui->actionBackgroundBlack)
        m_view->setBackgroundColor(Qt::black);
    else if(action == ui->actionBackgroundWhite)
        m_view->setBackgroundColor(Qt::white);
    else if(action == ui->actionBackgroundGray)
        m_view->setBackgroundColor(Qt::gray);
    else if(action == ui->actionBackgroundOther)
        m_view->setBackgroundColor(QColorDialog::getColor(m_view->backgroundColor(), this));

    m_view->update();
}

void BuilderWindow::setTool(QAction *action)
{
    if(action == ui->actionNavigate)
        setTool(navigateTool());
    else if(action == ui->actionBuild)
        setTool(buildTool());
    else if(action == ui->actionManipulate)
        setTool(manipulateTool());
}

void BuilderWindow::predictBonds()
{
    editor()->beginEdit();

    foreach(chemkit::Bond *bond, m_molecule->bonds())
        m_molecule->removeBond(bond);

    chemkit::BondPredictor predictor(m_molecule);

    QPair<chemkit::Atom *, chemkit::Atom *> bondedPair;
    foreach(bondedPair, predictor.predictedBonds()){
        editor()->addBond(bondedPair.first, bondedPair.second);
    }

    editor()->endEdit();
}

void BuilderWindow::adjustHydrogens()
{
    if(!m_molecule){
        return;
    }

    editor()->beginEdit();

    QSet<chemkit::Atom *> removedAtoms;

    foreach(chemkit::Atom *atom, m_molecule->atoms()){
        if(removedAtoms.contains(atom)){
            continue;
        }

        // add hydrogens
        while(atom->valence() < atom->expectedValence()){
            chemkit::Atom *hydrogen = editor()->addAtom(chemkit::Atom::Hydrogen);
            editor()->addBond(atom, hydrogen);
            editor()->setAtomPosition(hydrogen, atom->position().movedBy(chemkit::Vector::randomUnitVector()));
        }

        // remove hydrogens
        if(atom->valence() > atom->expectedValence()){
            foreach(chemkit::Atom *neighbor, atom->neighbors()){
                if(neighbor->isTerminalHydrogen()){
                    editor()->removeAtom(neighbor);
                    removedAtoms.insert(neighbor);

                    if(atom->valence() == atom->expectedValence()){
                        break;
                    }
                }
            }
        }
    }

    editor()->endEdit();
}

void BuilderWindow::setCanCut(bool canCut)
{
    ui->actionCut->setEnabled(canCut);
}

void BuilderWindow::setCanCopy(bool canCopy)
{
    ui->actionCopy->setEnabled(canCopy);
}

void BuilderWindow::setCanPaste(bool canPaste)
{
    ui->actionPaste->setEnabled(canPaste);
}

void BuilderWindow::setCanDelete(bool canDelete)
{
    ui->actionDelete->setEnabled(canDelete);
}

void BuilderWindow::moleculeProperties()
{
    MoleculePropertiesDialog dialog(m_molecule, this);
    dialog.exec();
}
