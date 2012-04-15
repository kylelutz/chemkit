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

#include "trajectoryviewerdemo.h"
#include "ui_trajectoryviewerdemo.h"

#include <chemkit/graphicsitem.h>
#include <chemkit/graphicscamera.h>
#include <chemkit/trajectoryfile.h>
#include <chemkit/graphicspainter.h>
#include <chemkit/trajectoryframe.h>
#include <chemkit/cartesiancoordinates.h>

// === GraphicsTrajectoryItem ============================================== //
class GraphicsTrajectoryItem : public chemkit::GraphicsItem
{
public:
    GraphicsTrajectoryItem();
    ~GraphicsTrajectoryItem();

    void setFrame(const chemkit::TrajectoryFrame *frame);
    const chemkit::TrajectoryFrame* frame() const;

    void paint(chemkit::GraphicsPainter *painter);

private:
    const chemkit::TrajectoryFrame *m_frame;
};

GraphicsTrajectoryItem::GraphicsTrajectoryItem()
{
    m_frame = 0;
}

GraphicsTrajectoryItem::~GraphicsTrajectoryItem()
{
}

void GraphicsTrajectoryItem::setFrame(const chemkit::TrajectoryFrame *frame)
{
    m_frame = frame;
}

const chemkit::TrajectoryFrame* GraphicsTrajectoryItem::frame() const
{
    return m_frame;
}

void GraphicsTrajectoryItem::paint(chemkit::GraphicsPainter *painter)
{
    if(!m_frame){
        return;
    }

    // set color to orange
    painter->setColor(QColor::fromRgb(255, 127, 0));

    // draw each atom as a sphere
    for(size_t i = 0; i < m_frame->size(); i++){
        const chemkit::Point3 position = m_frame->position(i);

        painter->drawSphere(position.cast<float>(), 0.1f);
    }
}

// === TrajectoryViewerDemo ================================================ //
TrajectoryViewerDemo::TrajectoryViewerDemo(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::TrajectoryViewerDemo)
{
    ui->setupUi(this);

    m_trajectoryItem = new GraphicsTrajectoryItem;
    ui->graphicsView->addItem(m_trajectoryItem);

    connect(ui->actionOpen, SIGNAL(triggered()), SLOT(openFile()));
    connect(ui->actionQuit, SIGNAL(triggered()), qApp, SLOT(quit()));
    connect(ui->openFileButton, SIGNAL(clicked()), ui->actionOpen, SLOT(trigger()));
    connect(ui->quitButton, SIGNAL(clicked()), ui->actionQuit, SLOT(trigger()));
    connect(ui->frameSlider, SIGNAL(valueChanged(int)), ui->frameSpinBox, SLOT(setValue(int)));
    connect(ui->frameSpinBox, SIGNAL(valueChanged(int)), ui->frameSlider, SLOT(setValue(int)));
    connect(ui->frameSpinBox, SIGNAL(valueChanged(int)), SLOT(setCurrentFrame(int)));
}

TrajectoryViewerDemo::~TrajectoryViewerDemo()
{
    delete ui;
}

void TrajectoryViewerDemo::setTrajectory(const boost::shared_ptr<chemkit::Trajectory> &trajectory)
{
    m_trajectory = trajectory;

    size_t frameCount = m_trajectory->frameCount();

    ui->frameSlider->setRange(1, frameCount);
    ui->frameSpinBox->setRange(1, frameCount);
    ui->frameCountLabel->setText(QString("/ %1").arg(frameCount));

    setCurrentFrame(1);

    if(m_trajectory && !m_trajectory->isEmpty()){
        const chemkit::TrajectoryFrame *frame = m_trajectory->frame(0);
        const chemkit::Point3 center = frame->coordinates()->center();

        ui->graphicsView->camera()->lookAt(center.cast<float>());
    }
}

void TrajectoryViewerDemo::openFile()
{
    QString fileName = QFileDialog::getOpenFileName(this,
                                                    tr("Open Trajectory File"),
                                                    0,
                                                    QString("All Files (*.*)"));

    if(!fileName.isEmpty()){
        openFile(fileName);
    }
}

void TrajectoryViewerDemo::openFile(const QString &fileName)
{
    std::string fileNameString = fileName.toStdString();

    chemkit::TrajectoryFile file(fileNameString);
    bool ok = file.read();
    if(!ok){
        QString errorString =
            QString("Failed to read file: %1").arg(file.errorString().c_str());
        QString("Failed to read file: %1").arg(file.errorString().c_str());

        QMessageBox::critical(this, "Read Error", errorString);
        return;
    }

    setTrajectory(file.trajectory());
}

void TrajectoryViewerDemo::setCurrentFrame(int index)
{
    if(!m_trajectory){
        return;
    }

    chemkit::TrajectoryFrame *frame = 0;

    if(static_cast<size_t>(index) - 1 < m_trajectory->frameCount()){
        frame = m_trajectory->frame(index - 1);
    }

    m_trajectoryItem->setFrame(frame);

    ui->graphicsView->update();
}

// === Main ================================================================ //
int main(int argc, char *argv[])
{
    QApplication app(argc, argv);

    TrajectoryViewerDemo window;
    window.show();

    if(app.argc() >= 2){
        window.openFile(app.arguments()[1]);
    }

    return app.exec();
}
