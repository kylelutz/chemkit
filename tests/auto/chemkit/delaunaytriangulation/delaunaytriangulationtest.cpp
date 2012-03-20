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

#include "delaunaytriangulationtest.h"

#include <algorithm>

#include <chemkit/delaunaytriangulation.h>

// This test case is based on the example presented on page 725 of the
// paper: "Three dimensional triganulations from local transformations"
// by Barry Joe (Siam J. Sci. Stat. Comput. Vol 10, No 4, 1989).
void DelaunayTriangulationTest::joe89()
{
    std::vector<chemkit::Point3> points;
    points.push_back(chemkit::Point3(0.054f, 0.099f, 0.993f));
    points.push_back(chemkit::Point3(0.066f, 0.756f, 0.910f));
    points.push_back(chemkit::Point3(0.076f, 0.578f, 0.408f));
    points.push_back(chemkit::Point3(0.081f, 0.036f, 0.954f));
    points.push_back(chemkit::Point3(0.082f, 0.600f, 0.726f));
    points.push_back(chemkit::Point3(0.085f, 0.327f, 0.731f));
    points.push_back(chemkit::Point3(0.123f, 0.666f, 0.842f));
    points.push_back(chemkit::Point3(0.161f, 0.303f, 0.975f));

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

    std::vector<std::vector<int> > tetrahedra = triangulation.tetrahedra();
    QCOMPARE(tetrahedra.size(), size_t(expectedTetrahedra.size()));

    for(unsigned int i = 0; i < tetrahedra.size(); i++){
        std::vector<int> tetrahedron = tetrahedra[i];
        std::sort(tetrahedron.begin(), tetrahedron.end());

        int foundCount = 0;

        for(int j = 0; j < expectedTetrahedra.size(); j++){
            QVector<int> expectedTetrahedron = expectedTetrahedra[j];

            if(tetrahedron == expectedTetrahedron.toStdVector()){
                foundCount++;
            }
        }

        if(foundCount != 1){
            qDebug() << "tetrahedron (index " << i << ") was not found. vertices: " << QVector<int>::fromStdVector(tetrahedron);
        }

        QCOMPARE(foundCount, 1);
    }
}

void DelaunayTriangulationTest::serine()
{
    // coordinates of the atoms
    std::vector<chemkit::Point3> points;
    points.push_back(chemkit::Point3(-0.1664, -1.0370, 0.4066));
    points.push_back(chemkit::Point3(1.2077, -0.5767, -0.0716));
    points.push_back(chemkit::Point3(-0.6079, -1.5894, -0.3173));
    points.push_back(chemkit::Point3(1.1440, -0.3456, -1.0571));
    points.push_back(chemkit::Point3(2.2495, -1.7077, 0.1008));
    points.push_back(chemkit::Point3(1.6659, 0.7153, 0.7175));
    points.push_back(chemkit::Point3(1.7844, 0.4727, 1.7759));
    points.push_back(chemkit::Point3(0.8959, 1.5129, 0.6034));
    points.push_back(chemkit::Point3(2.8918, 1.1700, 0.2007));
    points.push_back(chemkit::Point3(3.1444, 1.9558, 0.6711));
    points.push_back(chemkit::Point3(1.8101, -2.8570, 0.2804));
    points.push_back(chemkit::Point3(3.4579, -1.3878, 0.0035));
    points.push_back(chemkit::Point3(-0.0600, -1.6097, 1.2601));
    points.push_back(chemkit::Point3(-0.7527, -0.2118, 0.6162));

    // calculate delaunay triangulation
    chemkit::DelaunayTriangulation triangulation(points);
    QCOMPARE(triangulation.vertexCount(), 14);
    QCOMPARE(triangulation.edgeCount(), 60);
    QCOMPARE(triangulation.triangleCount(), 140);
    QCOMPARE(triangulation.tetrahedronCount(), 39);

    // weights (squared van der waals radii)
    std::vector<chemkit::Real> weights;
    weights.push_back(2.4025);
    weights.push_back(2.90);
    weights.push_back(1.44);
    weights.push_back(1.44);
    weights.push_back(2.90);
    weights.push_back(2.90);
    weights.push_back(1.44);
    weights.push_back(1.44);
    weights.push_back(2.3104);
    weights.push_back(1.44);
    weights.push_back(2.3104);
    weights.push_back(2.3104);
    weights.push_back(1.44);
    weights.push_back(1.44);

    // calculate weighted delaunay triangulation
    chemkit::DelaunayTriangulation weightedTriangulation(points, weights);
    QCOMPARE(weightedTriangulation.vertexCount(), 14);
    QCOMPARE(weightedTriangulation.edgeCount(), 60);
    QCOMPARE(weightedTriangulation.triangleCount(), 140);
    QCOMPARE(weightedTriangulation.tetrahedronCount(), 39);
}

QTEST_APPLESS_MAIN(DelaunayTriangulationTest)
