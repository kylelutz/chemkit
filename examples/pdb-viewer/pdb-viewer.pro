TEMPLATE = app
QT += opengl

INCLUDEPATH += ../../include
LIBS = -L../../lib -lchemkit -lchemkit-graphics

SOURCES += pdbviewerexample.cpp
HEADERS += pdbviewerexample.h
FORMS += pdbviewerexample.ui
