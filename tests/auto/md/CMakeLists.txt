if(NOT ${CHEMKIT_WITH_MD})
  return()
endif()

find_package(Chemkit COMPONENTS md)
include_directories(${CHEMKIT_INCLUDE_DIRS})

find_package(Qt4 4.6 COMPONENTS QtCore QtTest REQUIRED)
set(QT_DONT_USE_QTGUI TRUE)
set(QT_USE_QTTEST TRUE)
include(${QT_USE_FILE})

add_subdirectory(forcefield)
add_subdirectory(moleculegeometryoptimizer)
add_subdirectory(topology)
