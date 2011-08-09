find_package(Qt4 4.6 COMPONENTS QtCore REQUIRED)
set(QT_DONT_USE_QTGUI TRUE)
include(${QT_USE_FILE})

find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

include_directories(../../3rdparty/lemon)

set(SOURCES
  kekulizer.cpp
  smilesgraph.cpp
  smileslineformat.cpp
  smilesplugin.cpp
)

add_chemkit_plugin(smiles ${SOURCES})
target_link_libraries(smiles chemkit chemkit-io ${QT_LIBRARIES})
