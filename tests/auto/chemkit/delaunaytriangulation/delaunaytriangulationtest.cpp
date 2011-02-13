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

#include <QtTest>

#include <chemkit/delaunaytriangulation.h>

class DelaunayTriangulationTest : public QObject
{
    Q_OBJECT

    private slots:
        void joe89();
        void serine();
};

// This test case is based on the example presented on page 725 of the
// paper: "Three dimensional triganulations from local transformations"
// by Barry Joe (Siam J. Sci. Stat. Comput. Vol 10, No 4, 1989).
void DelaunayTriangulationTest::joe89()
{
    QVector<chemkit::Point> points;
    points.append(chemkit::Point(0.054f, 0.099f, 0.993f));
    points.append(chemkit::Point(0.066f, 0.756f, 0.910f));
    points.append(chemkit::Point(0.076f, 0.578f, 0.408f));
    points.append(chemkit::Point(0.081f, 0.036f, 0.954f));
    points.append(chemkit::Point(0.082f, 0.600f, 0.726f));
    points.append(chemkit::Point(0.085f, 0.327f, 0.731f));
    points.append(chemkit::Point(0.123f, 0.666f, 0.842f));
    points.append(chemkit::Point(0.161f, 0.303f, 0.975f));

    chemkit::DelaunayTriangulation triangulation(points);
    QCOMPARE(triangulation.vertexCount(), 8);
    QCOMPARE(triangulation.tetrahedronCount(), 13);

    QList<QVector<int> > expectedTetrahedra;
    expectedTetrahedra.append(QVector<int>() << 0 << 1 << 2 << 4);
    expectedTetrahedra.append(QVector<int>() << 0 << 1 << 4 << 5);
    expectedTetrahedra.append(QVector<int>() << 0 << 1 << 5 << 7);
    expectedTetrahedra.append(QVector<int>() << 0 << 2 << 3 << 5);
    expectedTetrahedra.append(QVector<int>() << 0 << 2 << 4 << 5);
    expectedTetrahedra.append(QVector<int>() << 0 << 3 << 5 << 7);
    expectedTetrahedra.append(QVector<int>() << 1 << 2 << 4 << 6);
    expectedTetrahedra.append(QVector<int>() << 1 << 4 << 5 << 7);
    expectedTetrahedra.append(QVector<int>() << 1 << 4 << 6 << 7);
    expectedTetrahedra.append(QVector<int>() << 2 << 3 << 5 << 7);
    expectedTetrahedra.append(QVector<int>() << 2 << 4 << 5 << 6);
    expectedTetrahedra.append(QVector<int>() << 2 << 5 << 6 << 7);
    expectedTetrahedra.append(QVector<int>() << 4 << 5 << 6 << 7);

    QList<QVector<int> > tetrahedra = triangulation.tetrahedra();
    QCOMPARE(tetrahedra.size(), expectedTetrahedra.size());

    for(int i = 0; i < tetrahedra.size(); i++){
        QVector<int> tetrahedron = tetrahedra[i];
        qSort(tetrahedron);

        int foundCount = 0;

        for(int j = 0; j < expectedTetrahedra.size(); j++){
            QVector<int> expectedTetrahedron = expectedTetrahedra[j];

            if(tetrahedron == expectedTetrahedron){
                foundCount++;
            }
        }

        if(foundCount != 1){
            qDebug() << "tetrahedron (index " << i << ") was not found. verticies: " << tetrahedron;
        }

        QCOMPARE(foundCount, 1);
    }
}

void DelaunayTriangulationTest::serine()
{
    // coordinates of the atoms
    QVector<chemkit::Point> points;
    points.append(chemkit::Point(-0.1664, -1.0370, 0.4066));
    points.append(chemkit::Point(1.2077, -0.5767, -0.0716));
    points.append(chemkit::Point(-0.6079, -1.5894, -0.3173));
    points.append(chemkit::Point(1.1440, -0.3456, -1.0571));
    points.append(chemkit::Point(2.2495, -1.7077, 0.1008));
    points.append(chemkit::Point(1.6659, 0.7153, 0.7175));
    points.append(chemkit::Point(1.7844, 0.4727, 1.7759));
    points.append(chemkit::Point(0.8959, 1.5129, 0.6034));
    points.append(chemkit::Point(2.8918, 1.1700, 0.2007));
    points.append(chemkit::Point(3.1444, 1.9558, 0.6711));
    points.append(chemkit::Point(1.8101, -2.8570, 0.2804));
    points.append(chemkit::Point(3.4579, -1.3878, 0.0035));
    points.append(chemkit::Point(-0.0600, -1.6097, 1.2601));
    points.append(chemkit::Point(-0.7527, -0.2118, 0.6162));

    // calculate delaunay triangulation
    chemkit::DelaunayTriangulation triangulation(points);
    QCOMPARE(triangulation.vertexCount(), 14);
    QCOMPARE(triangulation.edgeCount(), 60);
    QCOMPARE(triangulation.triangleCount(), 140);
    QCOMPARE(triangulation.tetrahedronCount(), 39);

    // weights (squared van der waals radii)
    QVector<chemkit::Float> weights;
    weights.append(2.4025);
    weights.append(2.90);
    weights.append(1.44);
    weights.append(1.44);
    weights.append(2.90);
    weights.append(2.90);
    weights.append(1.44);
    weights.append(1.44);
    weights.append(2.3104);
    weights.append(1.44);
    weights.append(2.3104);
    weights.append(2.3104);
    weights.append(1.44);
    weights.append(1.44);

    // calculate weighted delaunay triangulation
    chemkit::DelaunayTriangulation weightedTriangulation(points, weights);
    QCOMPARE(weightedTriangulation.vertexCount(), 14);
    QCOMPARE(weightedTriangulation.edgeCount(), 60);
    QCOMPARE(weightedTriangulation.triangleCount(), 140);
    QCOMPARE(weightedTriangulation.tetrahedronCount(), 39);
}

QTEST_APPLESS_MAIN(DelaunayTriangulationTest)
#include "delaunaytriangulationtest.moc"
