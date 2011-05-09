find_package(Qt4 4.6 COMPONENTS QtCore QtXml REQUIRED)
set(QT_DONT_USE_QTGUI TRUE)
set(QT_USE_QTXML TRUE)
include(${QT_USE_FILE})

find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

set(SOURCES
  pdbmlfileformat.cpp
  pdbmlplugin.cpp
)

set(MOC_HEADERS
  pdbmlplugin.h
)

qt4_wrap_cpp(MOC_SOURCES ${MOC_HEADERS})

add_chemkit_plugin(pdbml ${SOURCES} ${MOC_SOURCES})
target_link_libraries(pdbml chemkit ${QT_LIBRARIES})
