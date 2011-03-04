QT += testlib
CONFIG += console
INCLUDEPATH += ../../../../include/
LIBS = -L../../../../lib -lchemkit
SOURCES += atomtypertest.cpp \
    mockatomtyper.cpp
HEADERS += mockatomtyper.h
