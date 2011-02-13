include(../plugin.pri)

LIBS += -linchi

SOURCES += inchilineformat.cpp \
    inchikeylineformat.cpp \
    inchifileformat.cpp \
    inchiplugin.cpp
HEADERS += inchilineformat.h \
    inchikeylineformat.h \
    inchifileformat.h \
    inchiplugin.h
