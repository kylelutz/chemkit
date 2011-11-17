if(NOT ${CHEMKIT_WITH_IO})
  return()
endif()

find_package(Chemkit COMPONENTS io REQUIRED)
include_directories(${CHEMKIT_INCLUDE_DIRS})

find_package(Qt4 4.6 COMPONENTS QtCore QtXml REQUIRED)
set(QT_DONT_USE_QTGUI TRUE)
set(QT_USE_QTXML TRUE)
include(${QT_USE_FILE})

set(SOURCES
  cmlfileformat.cpp
  cmlplugin.cpp
)

add_chemkit_plugin(cml ${SOURCES})
target_link_libraries(cml ${CHEMKIT_LIBRARIES} ${QT_LIBRARIES})
