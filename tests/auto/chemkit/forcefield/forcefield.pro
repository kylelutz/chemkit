QT += testlib
CONFIG += console
INCLUDEPATH += ../../../../include/
LIBS = -L../../../../lib -lchemkit
SOURCES += forcefieldtest.cpp \
           mockforcefield.cpp
HEADERS += mockforcefield.h
