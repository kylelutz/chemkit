TEMPLATE = lib
TARGET = chemkit-widgets
DESTDIR = ../../lib
INCLUDEPATH += ../../include
LIBS += -L../../lib -lchemkit
DEFINES += CHEMKIT_WIDGETS_LIBRARY
SOURCES = moleculeeditor.cpp \
    periodictabledialog.cpp \
	periodictablewidget.cpp \
	widgets.cpp
HEADERS = moleculeeditor.h \
    periodictabledialog.h \
	periodictablewidget.h \
    widgets.h
FORMS = periodictabledialog.ui \
	periodictablewidget.ui
target.path = /usr/lib
INSTALLS += target
headers.files = $$HEADERS
headers.path = /usr/include/chemkit/
INSTALLS += headers
