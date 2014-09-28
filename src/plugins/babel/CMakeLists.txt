option(CHEMKIT_BUILD_PLUGIN_BABEL "Build the babel plugin." OFF)
if(NOT ${CHEMKIT_BUILD_PLUGIN_BABEL})
  return()
endif()

if(NOT ${CHEMKIT_WITH_IO})
  return()
endif()

find_package(Chemkit COMPONENTS io REQUIRED)
include_directories(${CHEMKIT_INCLUDE_DIRS})

find_package(Qt5 COMPONENTS Core REQUIRED)
set(QT_DONT_USE_QTGUI TRUE)
include_directories(${Qt5Core_INCLUDE_DIRS})

set(SOURCES
  babelfileformat.cpp
  babelplugin.cpp
)

add_chemkit_plugin(babel ${SOURCES})
target_link_libraries(babel Qt5::Core ${CHEMKIT_LIBRARIES})
