TEMPLATE = app
QT += opengl

INCLUDEPATH += ../../include
LIBS = -L../../lib -lchemkit -lchemkit-graphics

HEADERS += cubeviewerexample.h
SOURCES += cubeviewerexample.cpp
FORMS += cubeviewerexample.ui
