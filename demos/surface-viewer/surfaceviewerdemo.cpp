/******************************************************************************
**
** Copyright (C) 2009-2012 Kyle Lutz <kyle.r.lutz@gmail.com>
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

#include "surfaceviewerdemo.h"
#include "ui_surfaceviewerdemo.h"

#include <chemkit/moleculefile.h>
#include <chemkit/graphicscamera.h>
#include <chemkit/graphicsnavigationtool.h>

// === SurfaceViewerDemo =================================================== //
// --- Construction and Destruction ---------------------------------------- //
SurfaceViewerDemo::SurfaceViewerDemo(QWidget *parent)
    : QMainWindow(parent),
      ui(new Ui::SurfaceViewerDemo)
{
    // setup ui
    ui->setupUi(this);
    ui->colorModeComboBox->setCurrentIndex(4);

    // connect slots
    connect(ui->actionOpen, SIGNAL(triggered()), SLOT(openFile()));
    connect(ui->actionClose, SIGNAL(triggered()), SLOT(closeFile()));
    connect(ui->actionQuit, SIGNAL(triggered()), SLOT(quit()));
    connect(ui->surfaceTypeComboBox, SIGNAL(currentIndexChanged(int)), SLOT(surfaceTypeChanged(int)));
    connect(ui->colorModeComboBox, SIGNAL(currentIndexChanged(int)), SLOT(colorModeChanged(int)));
    connect(ui->opacitySlider, SIGNAL(valueChanged(int)), SLOT(opacitySliderChanged(int)));
    connect(ui->probeRadiusSpinBox, SIGNAL(valueChanged(double)), SLOT(probeRadiusChanged(double)));

    // setup graphics view
    m_moleculeItem = new chemkit::GraphicsMoleculeItem;
    ui->graphicsView->addItem(m_moleculeItem);

    m_molecularSurfaceItem = new chemkit::GraphicsMolecularSurfaceItem;
    m_molecularSurfaceItem->setSurfaceType(chemkit::MolecularSurface::VanDerWaals);
    ui->graphicsView->addItem(m_molecularSurfaceItem);
}

SurfaceViewerDemo::~SurfaceViewerDemo()
{
    closeFile();

    delete ui;
}

// --- Properties ---------------------------------------------------------- //
void SurfaceViewerDemo::setMolecule(const boost::shared_ptr<chemkit::Molecule> &molecule)
{
    m_moleculeItem->setMolecule(0);
    m_molecularSurfaceItem->setMolecule(0);

    m_molecule = molecule;

    if(molecule){
        m_moleculeItem->setMolecule(molecule.get());
        m_molecularSurfaceItem->setMolecule(molecule.get());
        ui->graphicsView->camera()->lookAt(molecule->center().cast<float>());
    }

    ui->graphicsView->update();
}

boost::shared_ptr<chemkit::Molecule> SurfaceViewerDemo::molecule() const
{
    return m_molecule;
}

// --- Slots --------------------------------------------------------------- //
void SurfaceViewerDemo::openFile()
{
    std::vector<std::string> formats = chemkit::MoleculeFile::formats();
    std::sort(formats.begin(), formats.end());

    QString formatsString;
    foreach(const std::string &format, formats){
        formatsString += QString("*.%1 ").arg(format.c_str());
    }

    QString fileName = QFileDialog::getOpenFileName(this,
                                                    "Open File",
                                                    0,
                                                    QString("Molecule Files (%1);;All Files (*.*)").arg(formatsString));

    if(!fileName.isEmpty()){
        openFile(fileName);
    }
}

void SurfaceViewerDemo::openFile(const QString &fileName)
{
    // close current file
    closeFile();

    // read file
    chemkit::MoleculeFile file(fileName.toStdString());
    if(!file.read()){
        QMessageBox::critical(this, "Error", QString("Error opening file: %1").arg(file.errorString().c_str()));
        return;
    }

    if(file.isEmpty()){
        QMessageBox::critical(this, "Error", QString("File is empty"));
        return;
    }

    setMolecule(file.molecule());
}

void SurfaceViewerDemo::closeFile()
{
    setMolecule(boost::shared_ptr<chemkit::Molecule>());
}

void SurfaceViewerDemo::quit()
{
    qApp->quit();
}

void SurfaceViewerDemo::surfaceTypeChanged(int index)
{
    if(index == 0){
        m_molecularSurfaceItem->setSurfaceType(chemkit::MolecularSurface::VanDerWaals);
    }
    else if(index == 1){
        m_molecularSurfaceItem->setSurfaceType(chemkit::MolecularSurface::SolventAccessible);
    }

    ui->graphicsView->update();
}

void SurfaceViewerDemo::colorModeChanged(int index)
{
    if(index == 0){
        m_molecularSurfaceItem->setColor(Qt::red);
        m_molecularSurfaceItem->setColorMode(chemkit::GraphicsMolecularSurfaceItem::SolidColor);
    }
    else if(index == 1){
        m_molecularSurfaceItem->setColor(Qt::green);
        m_molecularSurfaceItem->setColorMode(chemkit::GraphicsMolecularSurfaceItem::SolidColor);
    }
    else if(index == 2){
        m_molecularSurfaceItem->setColor(Qt::blue);
        m_molecularSurfaceItem->setColorMode(chemkit::GraphicsMolecularSurfaceItem::SolidColor);
    }
    else if(index == 3){
        m_molecularSurfaceItem->setColor(Qt::white);
        m_molecularSurfaceItem->setColorMode(chemkit::GraphicsMolecularSurfaceItem::SolidColor);
    }
    else if(index == 4){
        m_molecularSurfaceItem->setColorMode(chemkit::GraphicsMolecularSurfaceItem::AtomColor);
    }

    ui->graphicsView->update();
}

void SurfaceViewerDemo::opacitySliderChanged(int value)
{
    chemkit::Real opacity = value / 100.0;

    m_molecularSurfaceItem->setOpacity(opacity);

    ui->graphicsView->update();
}

void SurfaceViewerDemo::probeRadiusChanged(double radius)
{
    m_molecularSurfaceItem->setProbeRadius(radius);

    ui->graphicsView->update();
}

// === Main ================================================================ //
int main(int argc, char *argv[])
{
    QApplication app(argc, argv);

    SurfaceViewerDemo window;
    window.show();

    if(argc >= 2){
        window.openFile(argv[1]);
    }

    return app.exec();
}
