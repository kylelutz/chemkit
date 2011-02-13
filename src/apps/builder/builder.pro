TARGET = chemkit-builder
DESTDIR = ../../../bin
TEMPLATE = app
INCLUDEPATH += ../../../include
LIBS = -L../../../lib -lchemkit -lchemkit-graphics -lchemkit-widgets
QT += opengl
SOURCES += builder.cpp \
    buildertool.cpp \
    builderwindow.cpp \
    buildtool.cpp \
    displaysettingsdock.cpp \
    energyminimizationdock.cpp \
    energyminimizer.cpp \
    manipulatetool.cpp \
    moleculelistdock.cpp \
    moleculepropertiesdialog.cpp \
    navigatetool.cpp \
    toolsettingsdock.cpp
HEADERS += buildertool.h \
    builderwindow.h \
    buildtool.h \
    displaysettingsdock.h \
    energyminimizationdock.h \
    energyminimizer.h \
    manipulatetool.h \
    moleculelistdock.h \
    moleculepropertiesdialog.h \
    navigatetool.h \
    toolsettingsdock.h
FORMS += builderwindow.ui \
    displaysettingsdock.ui \
    energyminimizationdock.ui \
    moleculelistdock.ui \
    moleculepropertiesdialog.ui \
    toolsettingsdock.ui
target.path = /usr/bin
INSTALLS += target
