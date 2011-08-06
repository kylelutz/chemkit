find_package(Qt4 4.6 COMPONENTS QtCore REQUIRED)
set(QT_DONT_USE_QTGUI TRUE)
include(${QT_USE_FILE})

find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

set(SOURCES
  mopacplugin.cpp
  mopcrtfileformat.cpp
  mopinfileformat.cpp
)

add_chemkit_plugin(mopac ${SOURCES})
target_link_libraries(mopac chemkit ${QT_LIBRARIES})
