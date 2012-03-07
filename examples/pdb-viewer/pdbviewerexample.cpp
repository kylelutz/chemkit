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
        const boost::shared_ptr<chemkit::Polymer> &polymer = file->polymer();

        if(polymer){
            m_proteinItem->setPolymer(polymer.get());
            m_nucleicAcidItem->setPolymer(polymer.get());
            m_view->camera()->lookAt(polymer->center().cast<float>());
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

    QByteArray formatString = QFileInfo(fileName).suffix().toAscii();
    std::string format = formatString.constData();
    if(format == "xml")
        format = "pdbml";

    // open and read file
    chemkit::PolymerFile *file = new chemkit::PolymerFile;

    QByteArray fileNameString = fileName.toAscii();
    bool ok = file->read(fileNameString.constData(), format);
    if(!ok){
        QMessageBox::critical(this, "Error Reading File", file->errorString().c_str());
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
