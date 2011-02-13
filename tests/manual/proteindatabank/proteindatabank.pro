QT += testlib network
CONFIG += console
INCLUDEPATH += ../../../include/
LIBS = -L../../../lib -lchemkit -lchemkit-web
SOURCES = proteindatabanktest.cpp
