if(NOT ${CHEMKIT_WITH_WEB})
  return()
endif()

find_package(Chemkit COMPONENTS web REQUIRED)
include_directories(${CHEMKIT_INCLUDE_DIRS})

find_package(Qt4 4.6 COMPONENTS QtCore QtNetwork REQUIRED)
set(QT_DONT_USE_QTGUI TRUE)
include(${QT_USE_FILE})

set(SOURCES
  caslineformat.cpp
  casplugin.cpp
)

add_chemkit_plugin(cas ${SOURCES})
target_link_libraries(cas ${CHEMKIT_LIBRARIES} ${QT_LIBRARIES})
