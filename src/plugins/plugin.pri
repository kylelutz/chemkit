TEMPLATE = lib
CONFIG += plugin
QT -= gui

INCLUDEPATH += ../../../include
LIBS += -L../../../lib -lchemkit

DESTDIR = ../../../share/chemkit/plugins/

target.path = /usr/share/chemkit/plugins/
INSTALLS += target

# install target for plugin data files
exists($$PWD/$$TARGET/data) {
    pluginData.files = data/*
    pluginData.path = /usr/share/chemkit/plugins/data/$$TARGET/
    INSTALLS += pluginData

    # plugin data source directory
    PLUGIN_DATA_SRC = $$PWD/$$TARGET/data/

    # root plugin data directory
    PLUGIN_DATA_ROOT = $$DESTDIR/data/

    # destination directory for plugin data for each plugin target
    PLUGIN_DATA_DST = $$PLUGIN_DATA_ROOT/$$TARGET/

    # convert to proper path separators ('\' instead of '/') on windows
    win32-msvc* {
        PLUGIN_DATA_SRC = $$replace(PLUGIN_DATA_SRC, "/", "\\")
        PLUGIN_DATA_ROOT = $$replace(PLUGIN_DATA_ROOT, "/", "\\")
        PLUGIN_DATA_DST = $$replace(PLUGIN_DATA_DST, "/", "\\")
    }

    unix {
        COMMAND_SEPARATOR = ;
    }
    win32-msvc* {
        COMMAND_SEPARATOR = &
    }

    # create plugin root directory if it does not exist
    !exists($$PLUGIN_DATA_ROOT){
        QMAKE_POST_LINK += $$QMAKE_MKDIR $$PLUGIN_DATA_ROOT $$COMMAND_SEPARATOR
    }

    # create plugin data directory if it does not exist
    !exists($$PLUGIN_DATA_DST){
        QMAKE_POST_LINK += $$QMAKE_MKDIR $$PLUGIN_DATA_DST $$COMMAND_SEPARATOR
    }

    # copy plugin data files to plugin data directory
    QMAKE_POST_LINK += $$QMAKE_COPY $$PLUGIN_DATA_SRC* $$PLUGIN_DATA_DST
}
