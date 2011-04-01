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

#include "cubeviewerexample.h"
#include "ui_cubeviewerexample.h"

#include <chemkit/chemicalfile.h>
#include <chemkit/bondpredictor.h>
#include <chemkit/graphicsmaterial.h>
#include <chemkit/graphicsnavigationtool.h>

// === CubeViewerExample =================================================== //
// --- Construction and Destruction ---------------------------------------- //
CubeViewerExample::CubeViewerExample(QWidget *parent)
    : QMainWindow(parent),
      ui(new Ui::CubeViewerExample)
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

    // add navigation tool
    m_view->setTool(new chemkit::GraphicsNavigationTool);
}

CubeViewerExample::~CubeViewerExample()
{
    closeFile();

    delete ui;
}

// --- Slots --------------------------------------------------------------- //
void CubeViewerExample::openFile()
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

void CubeViewerExample::openFile(const QString &fileName)
{
    // close current file
    closeFile();

    // open new file
    chemkit::ChemicalFile file(fileName.toStdString());
    bool ok = file.read();
    if(!ok){
        QMessageBox::critical(this,
                              "Error Opening File",
                              QString("Failed to open file: %1").arg(file.errorString().c_str()));
        return;
    }

    // setup molecule
    chemkit::Molecule *molecule = file.molecule();
    chemkit::BondPredictor::predictBonds(molecule);
    m_moleculeItem->setMolecule(molecule);
    file.removeMolecule(molecule);

    // setup scalar fields and isosurface items
    m_positiveScalarField = readVolumeData(fileName);
    if(m_positiveScalarField){
        m_positiveSurfaceItem->setScalarField(m_positiveScalarField);
        m_positiveSurfaceItem->setPosition(m_positiveScalarField->origin());

        std::vector<chemkit::Float> values = m_positiveScalarField->data();
        for(unsigned int i = 0; i < values.size(); i++){
            values[i] = -values[i];
        }
        m_negativeScalarField = new chemkit::ScalarField(m_positiveScalarField->dimensions(),
                                                         m_positiveScalarField->cellDimensions(),
                                                         values);
        m_negativeSurfaceItem->setScalarField(m_negativeScalarField);
        m_negativeSurfaceItem->setPosition(m_positiveScalarField->origin());
    }

    // update the view
    m_view->update();
}

void CubeViewerExample::closeFile()
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

void CubeViewerExample::quit()
{
    qApp->quit();
}

void CubeViewerExample::isovalueChanged(int value)
{
    chemkit::Float isovalue = value / 1000.0;

    m_positiveSurfaceItem->setIsovalue(isovalue);
    m_negativeSurfaceItem->setIsovalue(isovalue);

    m_view->update();
}

void CubeViewerExample::opacityChanged(int value)
{
    chemkit::Float opacity = value / 100.0;

    m_positiveSurfaceItem->setOpacity(opacity);
    m_negativeSurfaceItem->setOpacity(opacity);

    m_view->update();
}

// --- Internal Methods ---------------------------------------------------- //
// Reads and returns the scalar field data from a .cube file.
chemkit::ScalarField* CubeViewerExample::readVolumeData(const QString &fileName) const
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
    std::vector<chemkit::Vector3f> axes(3);
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
    std::vector<chemkit::Float> volumeData;
    while(!file.atEnd()){
        line = file.readLine();
        lineItems = line.split(" ", QString::SkipEmptyParts);

        for(int i = 0; i < lineItems.size(); i++){
            volumeData.push_back(lineItems[i].toFloat());
        }
    }

    std::vector<chemkit::Float> cellLengths(3);
    for(int i = 0; i < 3; i++){
        cellLengths[i] = axes[i].length();
    }

    chemkit::ScalarField *scalarField = new chemkit::ScalarField(dimensions, cellLengths, volumeData);
    scalarField->setOrigin(origin);

    return scalarField;
}

// === Main ================================================================ //
int main(int argc, char *argv[])
{
    QApplication app(argc, argv);

    CubeViewerExample window;
    window.show();

    return app.exec();
}
