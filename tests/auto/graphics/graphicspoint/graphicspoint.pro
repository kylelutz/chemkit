QT += testlib opengl
CONFIG += console
INCLUDEPATH += ../../../../include/
LIBS = -L../../../../lib -lchemkit -lchemkit-graphics
SOURCES += graphicspointtest.cpp
