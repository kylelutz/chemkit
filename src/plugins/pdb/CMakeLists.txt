if(NOT ${CHEMKIT_WITH_IO})
  return()
endif()

find_package(Chemkit COMPONENTS io REQUIRED)
include_directories(${CHEMKIT_INCLUDE_DIRS})

find_package(Qt4 4.6 COMPONENTS QtCore REQUIRED)
set(QT_DONT_USE_QTGUI TRUE)
include(${QT_USE_FILE})

set(SOURCES
  pdbfileformat.cpp
  pdbplugin.cpp
)

add_chemkit_plugin(pdb ${SOURCES})
target_link_libraries(pdb ${CHEMKIT_LIBRARIES} ${QT_LIBRARIES})
