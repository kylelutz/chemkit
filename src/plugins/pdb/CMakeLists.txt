if(NOT ${CHEMKIT_WITH_IO})
  return()
endif()

find_package(Qt4 4.6 COMPONENTS QtCore REQUIRED)
set(QT_DONT_USE_QTGUI TRUE)
include(${QT_USE_FILE})

find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

set(SOURCES
  pdbfileformat.cpp
  pdbplugin.cpp
)

add_chemkit_plugin(pdb ${SOURCES})
target_link_libraries(pdb chemkit chemkit-io ${QT_LIBRARIES})
