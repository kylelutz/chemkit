find_package(Qt4 4.6 COMPONENTS QtCore QtGui QtTest REQUIRED)
set(QT_USE_QTTEST TRUE)
include(${QT_USE_FILE})

include_directories(../../../include)

add_subdirectory(moleculeeditor)
