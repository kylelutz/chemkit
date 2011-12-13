if(NOT ${CHEMKIT_WITH_WIDGETS})
  return()
endif()

find_package(Chemkit COMPONENTS widgets)
include_directories(${CHEMKIT_INCLUDE_DIRS})

find_package(Qt4 4.6 COMPONENTS QtCore QtGui QtTest REQUIRED)
set(QT_USE_QTTEST TRUE)
include(${QT_USE_FILE})

