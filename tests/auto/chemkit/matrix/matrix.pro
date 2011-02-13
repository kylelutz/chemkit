QT += testlib
CONFIG += console
INCLUDEPATH += ../../../../include/
LIBS = -L../../../../lib -lchemkit -lblas -llapack
SOURCES += matrixtest.cpp
