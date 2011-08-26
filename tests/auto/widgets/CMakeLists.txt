if(NOT ${CHEMKIT_WITH_WIDGETS})
  return()
endif()

find_package(Qt4 4.6 COMPONENTS QtCore QtGui QtTest REQUIRED)
set(QT_USE_QTTEST TRUE)
include(${QT_USE_FILE})

find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

include_directories(../../../include)

add_subdirectory(moleculeeditor)
