if(NOT ${CHEMKIT_WITH_WEB})
  return()
endif()

find_package(Qt4 4.6 COMPONENTS QtCore QtNetwork QtXml REQUIRED)
set(QT_DONT_USE_QTGUI TRUE)
set(QT_USE_QTNETWORK TRUE)
set(QT_USE_QTXML TRUE)
include(${QT_USE_FILE})

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

qt4_wrap_cpp(MOC_SOURCES ${MOC_HEADERS})

add_definitions(
  -DCHEMKIT_WEB_LIBRARY
)

add_chemkit_library(chemkit-web ${SOURCES} ${MOC_SOURCES})
target_link_libraries(chemkit-web ${CHEMKIT_LIBRARIES} ${QT_LIBRARIES})

# install header files
install(FILES ${HEADERS} DESTINATION include/chemkit/)
