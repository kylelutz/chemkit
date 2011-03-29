QT += testlib
CONFIG += console
INCLUDEPATH += ../../../../include/
LIBS = -L../../../../lib -lchemkit
HEADERS = mockclass.h \
    mockplugin.h
SOURCES = mockplugin.cpp \
    plugintest.cpp

