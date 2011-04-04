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

#include "plugintest.h"

#include <string>
#include <vector>
#include <algorithm>

#include <chemkit/plugin.h>
#include <chemkit/pluginmanager.h>

#include "mockclass.h"
#include "mockplugin.h"

void PluginTest::initTestCase()
{
    m_plugin = new MockPlugin();
}

void PluginTest::name()
{
    QCOMPARE(m_plugin->name(), std::string("mock"));
}

void PluginTest::registerClass()
{
    bool ok = m_plugin->registerClass("mockplugin");
    QVERIFY(ok);

    std::vector<std::string> plugins = chemkit::PluginManager::instance()->pluginClassNames<MockClass>();
    QVERIFY(plugins.size() == 1);
    QVERIFY(std::find(plugins.begin(), plugins.end(), "mockplugin") != plugins.end());

    ok = m_plugin->unregisterClass("mockplugin");
    QVERIFY(ok);

    plugins = chemkit::PluginManager::instance()->pluginClassNames<MockClass>();
    QVERIFY(plugins.size() == 0);
    QVERIFY(std::find(plugins.begin(), plugins.end(), "mockplugin") == plugins.end());
}

void PluginTest::cleanupTestCase()
{
    delete m_plugin;
}

QTEST_APPLESS_MAIN(PluginTest)
