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

#include "graphicsscenetest.h"

#include <chemkit/graphicsitem.h>
#include <chemkit/graphicsview.h>
#include <chemkit/graphicsscene.h>

void GraphicsSceneTest::items()
{
    chemkit::GraphicsScene *scene = new chemkit::GraphicsScene;
    QCOMPARE(scene->itemCount(), 0);
    QCOMPARE(scene->size(), 0);
    QCOMPARE(scene->isEmpty(), true);
    QCOMPARE(scene->items().size(), 0);

    chemkit::GraphicsItem *item = new chemkit::GraphicsItem(0);
    scene->addItem(item);
    QCOMPARE(scene->itemCount(), 1);
    QCOMPARE(scene->isEmpty(), false);
    QVERIFY(scene->items()[0] == item);

    chemkit::GraphicsItem *item2 = new chemkit::GraphicsItem(0);
    scene->addItem(item2);
    QCOMPARE(scene->itemCount(), 2);
    QVERIFY(scene->items()[1] == item2);

    scene->deleteItem(item);
    QCOMPARE(scene->itemCount(), 1);
    QVERIFY(scene->item(0) == item2);

    scene->deleteItem(item2);
    QCOMPARE(scene->itemCount(), 0);
    QCOMPARE(scene->isEmpty(), true);
}

void GraphicsSceneTest::views()
{
    chemkit::GraphicsScene *scene = new chemkit::GraphicsScene;
    QCOMPARE(scene->views().size(), 0);

    chemkit::GraphicsView *view = new chemkit::GraphicsView(scene);
    QCOMPARE(scene->views().size(), 1);
    QVERIFY(scene->views()[0] == view);

    chemkit::GraphicsView *view2 = new chemkit::GraphicsView(scene);
    QCOMPARE(scene->views().size(), 2);
    QVERIFY(scene->views()[1] == view2);

    delete view;
    QCOMPARE(scene->views().size(), 1);
    QVERIFY(scene->views()[0] == view2);

    delete view2;
    QCOMPARE(scene->views().size(), 0);

    delete scene;
}

QTEST_MAIN(GraphicsSceneTest)
