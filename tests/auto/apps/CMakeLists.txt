if(NOT ${CHEMKIT_BUILD_APPS})
  return()
endif()

find_package(Qt4 4.6 COMPONENTS QtCore QtTest REQUIRED)
set(QT_DONT_USE_QTGUI TRUE)
set(QT_USE_QTTEST TRUE)
include(${QT_USE_FILE})

add_subdirectory(convert)
add_subdirectory(grep)
add_subdirectory(gen3d)
add_subdirectory(translate)
