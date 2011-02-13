QT += testlib
CONFIG += console
INCLUDEPATH += ../../../../include/
LIBS = -L../../../../lib \
    -lchemkit
SOURCES += moleculardescriptortest.cpp \
    mockdescriptor.cpp
HEADERS += mockdescriptor.h
