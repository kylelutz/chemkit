TARGET = chemkit-convert
DESTDIR = ../../../bin

TEMPLATE = app
CONFIG += console
INCLUDEPATH += ../../../include
LIBS = -L../../../lib -lchemkit

SOURCES += convert.cpp

target.path = /usr/bin
INSTALLS += target
