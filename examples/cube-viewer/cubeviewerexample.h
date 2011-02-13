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
        chemkit::GraphicsView *m_view;
        chemkit::GraphicsMoleculeItem *m_moleculeItem;
        chemkit::GraphicsIsosurfaceItem *m_positiveSurfaceItem;
        chemkit::GraphicsIsosurfaceItem *m_negativeSurfaceItem;
        chemkit::ScalarField *m_positiveScalarField;
        chemkit::ScalarField *m_negativeScalarField;
};

#endif // CUBEVIEWEREXAMPLE_H
