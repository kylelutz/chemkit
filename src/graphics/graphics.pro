TEMPLATE = lib
TARGET = chemkit-graphics
DESTDIR = ../../lib
INCLUDEPATH += ../../include
QT += opengl
LIBS += -L../../lib \
    -lchemkit \
    -lblas \
    -llapack
DEFINES += CHEMKIT_GRAPHICS_LIBRARY
SOURCES += graphics.cpp \
    graphicsatomcolormap.cpp \
    graphicsatomitem.cpp \
    graphicsbonditem.cpp \
    graphicscamera.cpp \
    graphicscylinder.cpp \
    graphicscylinderitem.cpp \
    graphicsisosurfaceitem.cpp \
    graphicsitem.cpp \
    graphicslight.cpp \
    graphicsmaterial.cpp \
    graphicsmolecularsurfaceitem.cpp \
    graphicsmoleculeitem.cpp \
    graphicsnavigationtool.cpp \
    graphicsnucleicaciditem.cpp \
    graphicsoverlay.cpp \
    graphicspainter.cpp \
    graphicsproteincoilitem.cpp \
    graphicsproteinhelixitem.cpp \
    graphicsproteinitem.cpp \
    graphicsproteinsheetitem.cpp \
    graphicsray.cpp \
    graphicsringitem.cpp \
    graphicsscene.cpp \
    graphicssphere.cpp \
    graphicssphereitem.cpp \
    graphicstool.cpp \
    graphicstransform.cpp \
    graphicsview.cpp \
    graphicsvertexbuffer.cpp
HEADERS += graphics.h \
    graphicsatomcolormap.h \
    graphicsatomitem.h \
    graphicsbonditem.h \
    graphicscamera.h \
    graphicscylinder.h \
    graphicscylinderitem.h \
    graphicsisosurfaceitem.h \
    graphicsitem.h \
    graphicslight.h \
    graphicsmaterial.h \
    graphicsmolecularsurfaceitem.h \
    graphicsmoleculeitem.h \
    graphicsnavigationtool.h \
    graphicsnucleicaciditem.h \
    graphicsoverlay.h \
    graphicspainter.h \
    graphicsproteincoilitem.h \
    graphicsproteinhelixitem.h \
    graphicsproteinitem.h \
    graphicsproteinsheetitem.h \
    graphicsray.h \
    graphicsringitem.h \
    graphicsscene.h \
    graphicssphere.h \
    graphicssphereitem.h \
    graphicstool.h \
    graphicstransform.h \
    graphicsview.h \
    graphicsvertexbuffer.h
RESOURCES = shaders.qrc
OTHER_FILES += shaders/phong.frag \
    shaders/phong.vert \
    shaders/flat.frag \
    shaders/flat.vert
target.path = /usr/lib
INSTALLS += target
headers.files = $$HEADERS
headers.path = /usr/include/chemkit/
INSTALLS += headers
