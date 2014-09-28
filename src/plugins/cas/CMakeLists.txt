if(NOT ${CHEMKIT_WITH_WEB})
  return()
endif()

find_package(Chemkit COMPONENTS web REQUIRED)
include_directories(${CHEMKIT_INCLUDE_DIRS})

find_package(Qt5 COMPONENTS Core Network REQUIRED)
set(QT_DONT_USE_QTGUI TRUE)
include_directories(${Qt5Core_INCLUDE_DIRS})
include_directories(${Qt5Network_INCLUDE_DIRS})

set(SOURCES
  caslineformat.cpp
  casplugin.cpp
)

add_chemkit_plugin(cas ${SOURCES})
target_link_libraries(cas Qt5::Core Qt5::Network ${CHEMKIT_LIBRARIES})
