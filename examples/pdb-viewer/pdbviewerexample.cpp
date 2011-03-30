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

#include "pdbviewerexample.h"
#include "ui_pdbviewerexample.h"

#include <chemkit/polymer.h>
#include <chemkit/bondpredictor.h>
#include <chemkit/graphicscamera.h>
#include <chemkit/graphicsnavigationtool.h>

// === PdbViewerWindow ===================================================== //
// --- Construction and Destruction ---------------------------------------- //
PdbViewerWindow::PdbViewerWindow(QWidget *parent)
    : QMainWindow(parent),
      ui(new Ui::PdbViewerWindow)
{
    m_file = 0;

    // setup ui
    ui->setupUi(this);

    // setup icons for the menu actions
    ui->actionOpen->setIcon(QIcon::fromTheme("document-open", style()->standardIcon(QStyle::SP_DialogOpenButton)));
    ui->actionQuit->setIcon(QIcon::fromTheme("application-exit", style()->standardIcon(QStyle::SP_DialogCloseButton)));

    connect(ui->actionOpen, SIGNAL(triggered()), SLOT(openFile()));
    connect(ui->actionQuit, SIGNAL(triggered()), SLOT(quit()));

    // setup view
    m_view = new chemkit::GraphicsView;
    setCentralWidget(m_view);
    m_view->setTool(new chemkit::GraphicsNavigationTool);

    m_proteinItem = new chemkit::GraphicsProteinItem;
    m_view->addItem(m_proteinItem);
    m_nucleicAcidItem = new chemkit::GraphicsNucleicAcidItem;
    m_view->addItem(m_nucleicAcidItem);
}

PdbViewerWindow::~PdbViewerWindow()
{
    setFile(0);

    delete m_view;
}

// --- Properties ---------------------------------------------------------- //
void PdbViewerWindow::setFile(chemkit::PolymerFile *file)
{
    if(m_file){
        delete m_file;
        m_proteinItem->setPolymer(0);
        m_nucleicAcidItem->setPolymer(0);
    }

    m_file = file;

    if(file){
        chemkit::Polymer *polymer = file->polymer();

        if(polymer){
            m_proteinItem->setPolymer(polymer);
            m_nucleicAcidItem->setPolymer(polymer);
            m_view->camera()->lookAt(polymer->center());
        }
    }

    m_view->update();
}

// --- Slots --------------------------------------------------------------- //
void PdbViewerWindow::openFile()
{
    QString fileName = QFileDialog::getOpenFileName(this,
                                                    tr("Open File"),
                                                    0,
                                                    QString("PDB File (*.pdb);;All Files (*.*)"));

    if(!fileName.isEmpty()){
        openFile(fileName);
    }
}

void PdbViewerWindow::openFile(const QString &fileName)
{
    // close current file
    setFile(0);

    std::string format = QFileInfo(fileName).suffix().toStdString();
    if(format == "xml")
        format = "pdbml";

    // open and read file
    chemkit::PolymerFile *file = new chemkit::PolymerFile;

    bool ok = file->read(fileName, format);
    if(!ok){
        QMessageBox::critical(this, "Error Reading File", file->errorString());
        delete file;
        return;
    }

    setFile(file);
}

void PdbViewerWindow::quit()
{
    qApp->exit();
}

// === Main ================================================================ //
int main(int argc, char *argv[])
{
    QApplication app(argc, argv);

    PdbViewerWindow window;
    window.show();

    if(app.argc() >= 2)
        window.openFile(app.arguments()[1]);

    return app.exec();
}
