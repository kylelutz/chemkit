if(NOT ${CHEMKIT_WITH_WEB})
  return()
endif()

find_package(Qt5 COMPONENTS Core Network Xml REQUIRED)
set(QT_DONT_USE_QTGUI TRUE)
set(QT_USE_QTNETWORK TRUE)
set(QT_USE_QTXML TRUE)
include_directories(${Qt5Core_INCLUDE_DIRS})
include_directories(${Qt5Network_INCLUDE_DIRS})
include_directories(${Qt5Xml_INCLUDE_DIRS})

find_package(Chemkit COMPONENTS io REQUIRED)
include_directories(${CHEMKIT_INCLUDE_DIRS})

set(HEADERS
  proteindatabank.h
  pubchem.h
  web.h
)

set(SOURCES
  downloadthread.cpp
  proteindatabank.cpp
  pubchem.cpp
  pubchemquery.cpp
  pubchemquerythread.cpp
  web.cpp
)

set(MOC_HEADERS
  downloadthread.h
  pubchemquerythread.h
)

qt5_wrap_cpp(MOC_SOURCES ${MOC_HEADERS})

add_definitions(
  -DCHEMKIT_WEB_LIBRARY
)

add_chemkit_library(chemkit-web ${SOURCES} ${MOC_SOURCES})
target_link_libraries(chemkit-web Qt5::Core Qt5::Network Qt5::Xml ${CHEMKIT_LIBRARIES})

# install header files
install(FILES ${HEADERS} DESTINATION include/chemkit/)
