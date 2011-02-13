TEMPLATE = lib
TARGET = chemkit-web
DESTDIR = ../../lib
INCLUDEPATH += ../../include
LIBS = -L../../lib -lchemkit
QT += network \
    xml
QT -= gui
DEFINES += CHEMKIT_WEB_LIBRARY
SOURCES = downloadthread.cpp \
    proteindatabank.cpp \
    pubchem.cpp \
    pubchemquery.cpp \
    pubchemquerythread.cpp \
    web.cpp
HEADERS = downloadthread.h \
    proteindatabank.h \
    pubchem.h \
    pubchemquery.h \
    pubchemquerythread.h \
    web.h
target.path = /usr/lib
INSTALLS += target
headers.files = proteindatabank.h \
    pubchem.h \
    web.h
headers.path = /usr/include/chemkit/
INSTALLS += headers
