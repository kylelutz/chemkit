if(NOT ${CHEMKIT_WITH_IO})
  return()
endif()

find_package(Chemkit COMPONENTS io REQUIRED)
include_directories(${CHEMKIT_INCLUDE_DIRS})

find_package(Qt4 4.6 COMPONENTS QtCore REQUIRED)
set(QT_DONT_USE_QTGUI TRUE)
include(${QT_USE_FILE})

set(SOURCES
  fhzfileformat.cpp
  fhzplugin.cpp
)

add_chemkit_plugin(fhz ${SOURCES})
target_link_libraries(fhz ${CHEMKIT_LIBRARIES} ${QT_LIBRARIES})
