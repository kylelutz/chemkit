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

#include "cubeviewerdemo.h"
#include "ui_cubeviewerdemo.h"

#include <chemkit/moleculefile.h>
#include <chemkit/bondpredictor.h>
#include <chemkit/graphicsmaterial.h>
#include <chemkit/graphicsnavigationtool.h>

// === CubeViewerDemo =================================================== //
// --- Construction and Destruction ---------------------------------------- //
CubeViewerDemo::CubeViewerDemo(QWidget *parent)
    : QMainWindow(parent),
      ui(new Ui::CubeViewerDemo)
{
    // setup ui
    ui->setupUi(this);

    // setup graphics widget widget
    m_view = new chemkit::GraphicsView;
    setCentralWidget(m_view);

    // connect signals to slots
    connect(ui->actionOpen, SIGNAL(triggered()), SLOT(openFile()));
    connect(ui->actionClose, SIGNAL(triggered()), SLOT(closeFile()));
    connect(ui->actionQuit, SIGNAL(triggered()), SLOT(quit()));
    connect(ui->isovalueSlider, SIGNAL(valueChanged(int)), SLOT(isovalueChanged(int)));
    connect(ui->opacitySlider, SIGNAL(valueChanged(int)), SLOT(opacityChanged(int)));

    // setup graphics items
    m_moleculeItem = new chemkit::GraphicsMoleculeItem;
    m_view->addItem(m_moleculeItem);

    m_positiveSurfaceItem = new chemkit::GraphicsIsosurfaceItem;
    m_positiveSurfaceItem->setColor(Qt::red);
    m_positiveSurfaceItem->material()->setSpecularColor(Qt::transparent);
    m_view->addItem(m_positiveSurfaceItem);

    m_negativeSurfaceItem = new chemkit::GraphicsIsosurfaceItem;
    m_negativeSurfaceItem->setColor(Qt::blue);
    m_negativeSurfaceItem->material()->setSpecularColor(Qt::transparent);
    m_view->addItem(m_negativeSurfaceItem);

    // setup scalar fields
    m_positiveScalarField = 0;
    m_negativeScalarField = 0;
}

CubeViewerDemo::~CubeViewerDemo()
{
    closeFile();

    delete ui;
}

// --- Slots --------------------------------------------------------------- //
void CubeViewerDemo::openFile()
{
    // close current file
    closeFile();

    QString fileName = QFileDialog::getOpenFileName(this,
                                                    "Open File",
                                                    QString(),
                                                    QString("Cube Files (*.cube)"));
    if(!fileName.isEmpty()){
        openFile(fileName);
    }
}

void CubeViewerDemo::openFile(const QString &fileName)
{
    // close current file
    closeFile();

    // open new file
    QByteArray fileNameString = fileName.toAscii();
    chemkit::MoleculeFile file(fileNameString.constData());
    bool ok = file.read();
    if(!ok){
        QMessageBox::critical(this,
                              "Error Opening File",
                              QString("Failed to open file: %1").arg(file.errorString().c_str()));
        return;
    }

    // setup molecule
    m_molecule = file.molecule();
    chemkit::BondPredictor::predictBonds(m_molecule.get());
    m_moleculeItem->setMolecule(m_molecule.get());

    // setup scalar fields and isosurface items
    m_positiveScalarField = readVolumeData(fileName);
    if(m_positiveScalarField){
        m_positiveSurfaceItem->setScalarField(m_positiveScalarField);
        m_positiveSurfaceItem->setPosition(m_positiveScalarField->origin().cast<float>());

        std::vector<chemkit::Real> values = m_positiveScalarField->data();
        for(unsigned int i = 0; i < values.size(); i++){
            values[i] = -values[i];
        }
        m_negativeScalarField = new chemkit::ScalarField(m_positiveScalarField->dimensions(),
                                                         m_positiveScalarField->cellDimensions(),
                                                         values);
        m_negativeSurfaceItem->setScalarField(m_negativeScalarField);
        m_negativeSurfaceItem->setPosition(m_positiveScalarField->origin().cast<float>());
    }

    // update the view
    m_view->update();
}

void CubeViewerDemo::closeFile()
{
    if(m_positiveScalarField){
        delete m_positiveScalarField;
        m_positiveScalarField = 0;
    }

    if(m_negativeScalarField){
        delete m_negativeScalarField;
        m_negativeScalarField = 0;
    }
}

void CubeViewerDemo::quit()
{
    qApp->quit();
}

void CubeViewerDemo::isovalueChanged(int value)
{
    chemkit::Real isovalue = value / 1000.0;

    m_positiveSurfaceItem->setIsovalue(isovalue);
    m_negativeSurfaceItem->setIsovalue(isovalue);

    m_view->update();
}

void CubeViewerDemo::opacityChanged(int value)
{
    chemkit::Real opacity = value / 100.0;

    m_positiveSurfaceItem->setOpacity(opacity);
    m_negativeSurfaceItem->setOpacity(opacity);

    m_view->update();
}

// --- Internal Methods ---------------------------------------------------- //
// Reads and returns the scalar field data from a .cube file.
chemkit::ScalarField* CubeViewerDemo::readVolumeData(const QString &fileName) const
{
    QFile file(fileName);
    bool ok = file.open(QFile::ReadOnly);
    if(!ok){
        qDebug() << "Error: failed to read cube file: " << file.errorString();
        return 0;
    }

    // title line
    QString line = file.readLine();

    // comment line
    line = file.readLine();

    // atom count and origin coordinates line
    line = file.readLine();
    QStringList lineItems = line.split(" ", QString::SkipEmptyParts);
    if(lineItems.size() < 4){
        qDebug() << "Error: Cube file counts line too short.";
        return 0;
    }

    bool negativeAtomCount = false;
    int atomCount = lineItems[0].toInt();
    if(atomCount < 0){
        negativeAtomCount = true;
        atomCount = qAbs(atomCount);
    }

    chemkit::Point3 origin(lineItems[1].toDouble(),
                           lineItems[2].toDouble(),
                           lineItems[3].toDouble());

    // voxel count and axes
    std::vector<int> dimensions(3);
    std::vector<chemkit::Vector3> axes(3);
    for(int i = 0; i < 3; i++){
        line = file.readLine();
        lineItems = line.split(" ", QString::SkipEmptyParts);

        if(lineItems.size() < 4){
            continue;
        }

        dimensions[i] = lineItems[0].toInt();
        axes[i] = chemkit::Vector3(lineItems[1].toDouble(),
                                   lineItems[2].toDouble(),
                                   lineItems[3].toDouble());
    }

    // read past atoms
    for(int i = 0; i < atomCount; i++){
        line = file.readLine();
    }

    // a negative atom count indicates that the next line will
    // contain the orbital count and orbital number
    int orbitalCount = 0;
    int orbitalNumber = 0;
    if(negativeAtomCount){
        line = file.readLine();
        lineItems = line.split(" ", QString::SkipEmptyParts);
        if(lineItems.size() >= 2){
            orbitalCount = lineItems[0].toInt();
            orbitalNumber = lineItems[1].toInt();
        }
    }

    // read volume data
    std::vector<chemkit::Real> volumeData;
    while(!file.atEnd()){
        line = file.readLine();
        lineItems = line.split(" ", QString::SkipEmptyParts);

        for(int i = 0; i < lineItems.size(); i++){
            volumeData.push_back(lineItems[i].toFloat());
        }
    }

    std::vector<chemkit::Real> cellLengths(3);
    for(int i = 0; i < 3; i++){
        cellLengths[i] = axes[i].norm();
    }

    chemkit::ScalarField *scalarField = new chemkit::ScalarField(dimensions, cellLengths, volumeData);
    scalarField->setOrigin(origin);

    return scalarField;
}

// === Main ================================================================ //
int main(int argc, char *argv[])
{
    QApplication app(argc, argv);

    CubeViewerDemo window;
    window.show();

    return app.exec();
}
