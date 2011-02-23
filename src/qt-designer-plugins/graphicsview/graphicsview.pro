TEMPLATE = lib
CONFIG += designer plugin
QT += opengl
DESTDIR = ../../../lib/qt-designer-plugins
INCLUDEPATH += ../../../include
LIBS = -L../../../lib -lchemkit-graphics

SOURCES = graphicsviewdesignerplugin.cpp
HEADERS = graphicsviewdesignerplugin.h

target.path = $$[QT_INSTALL_PLUGINS]/designer
INSTALLS += target
