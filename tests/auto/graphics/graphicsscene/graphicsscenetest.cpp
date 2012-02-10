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

#include "graphicsscenetest.h"

#include <chemkit/graphicsitem.h>
#include <chemkit/graphicsview.h>
#include <chemkit/graphicsscene.h>

void GraphicsSceneTest::items()
{
    boost::shared_ptr<chemkit::GraphicsScene> scene(new chemkit::GraphicsScene);
    QCOMPARE(scene->itemCount(), size_t(0));
    QCOMPARE(scene->size(), size_t(0));
    QCOMPARE(scene->isEmpty(), true);
    QCOMPARE(scene->items().size(), size_t(0));

    chemkit::GraphicsItem *item = new chemkit::GraphicsItem(0);
    scene->addItem(item);
    QCOMPARE(scene->itemCount(), size_t(1));
    QCOMPARE(scene->isEmpty(), false);
    QVERIFY(scene->items()[0] == item);

    chemkit::GraphicsItem *item2 = new chemkit::GraphicsItem(0);
    scene->addItem(item2);
    QCOMPARE(scene->itemCount(), size_t(2));
    QVERIFY(scene->items()[1] == item2);

    scene->deleteItem(item);
    QCOMPARE(scene->itemCount(), size_t(1));
    QVERIFY(scene->item(0) == item2);

    scene->deleteItem(item2);
    QCOMPARE(scene->itemCount(), size_t(0));
    QCOMPARE(scene->isEmpty(), true);
}

void GraphicsSceneTest::views()
{
    boost::shared_ptr<chemkit::GraphicsScene> scene(new chemkit::GraphicsScene);
    QCOMPARE(scene->views().size(), size_t(0));

    chemkit::GraphicsView *view = new chemkit::GraphicsView(scene);
    QCOMPARE(scene->views().size(), size_t(1));
    QVERIFY(scene->views()[0] == view);

    chemkit::GraphicsView *view2 = new chemkit::GraphicsView(scene);
    QCOMPARE(scene->views().size(), size_t(2));
    QVERIFY(scene->views()[1] == view2);

    delete view;
    QCOMPARE(scene->views().size(), size_t(1));
    QVERIFY(scene->views()[0] == view2);

    delete view2;
    QCOMPARE(scene->views().size(), size_t(0));
}

QTEST_MAIN(GraphicsSceneTest)
