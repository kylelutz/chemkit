TARGET = chemkit-grep
DESTDIR = ../../../bin

TEMPLATE = app
CONFIG += console
INCLUDEPATH += ../../../include
LIBS = -L../../../lib -lchemkit 

SOURCES += grep.cpp

target.path = /usr/bin
INSTALLS += target
