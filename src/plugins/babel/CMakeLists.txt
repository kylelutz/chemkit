option(CHEMKIT_BUILD_PLUGIN_BABEL "Build the babel plugin." OFF)
if(NOT ${CHEMKIT_BUILD_PLUGIN_BABEL})
  return()
endif()

if(NOT ${CHEMKIT_WITH_IO})
  return()
endif()

find_package(Chemkit COMPONENTS io REQUIRED)
include_directories(${CHEMKIT_INCLUDE_DIRS})

find_package(Qt4 4.6 COMPONENTS QtCore REQUIRED)
set(QT_DONT_USE_QTGUI TRUE)
include(${QT_USE_FILE})

set(SOURCES
  babelfileformat.cpp
  babelplugin.cpp
)

add_chemkit_plugin(babel ${SOURCES})
target_link_libraries(babel ${CHEMKIT_LIBRARIES} ${QT_LIBRARIES})
