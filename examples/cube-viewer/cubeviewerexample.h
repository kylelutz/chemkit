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

#ifndef CUBEVIEWEREXAMPLE_H
#define CUBEVIEWEREXAMPLE_H

#include <QtGui>

#include <chemkit/scalarfield.h>
#include <chemkit/graphicsview.h>
#include <chemkit/graphicsmoleculeitem.h>
#include <chemkit/graphicsisosurfaceitem.h>

namespace Ui {
    class CubeViewerExample;
}

class CubeViewerExample : public QMainWindow
{
    Q_OBJECT

    public:
        CubeViewerExample(QWidget *parent = 0);
        ~CubeViewerExample();

    public slots:
        void openFile();
        void openFile(const QString &fileName);
        void closeFile();
        void quit();
        void opacityChanged(int value);
        void isovalueChanged(int value);

    private:
        chemkit::ScalarField* readVolumeData(const QString &fileName) const;

    private:
        Ui::CubeViewerExample *ui;
        boost::shared_ptr<chemkit::Molecule> m_molecule;
        chemkit::GraphicsView *m_view;
        chemkit::GraphicsMoleculeItem *m_moleculeItem;
        chemkit::GraphicsIsosurfaceItem *m_positiveSurfaceItem;
        chemkit::GraphicsIsosurfaceItem *m_negativeSurfaceItem;
        chemkit::ScalarField *m_positiveScalarField;
        chemkit::ScalarField *m_negativeScalarField;
};

#endif // CUBEVIEWEREXAMPLE_H
