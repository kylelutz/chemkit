TARGET = inchi
TEMPLATE = lib
DESTDIR = ../../../lib
CONFIG += staticlib
QT -= core gui

# The inchi sources are filled with warnings. Turn
# them off get rid of the noise.
CONFIG += warn_off

SOURCES = *.c
HEADERS = *.h
