find_package(Qt4 4.6 COMPONENTS QtCore REQUIRED)
set(QT_DONT_USE_QTGUI TRUE)
include(${QT_USE_FILE})

find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

set(SOURCES
  xyzfileformat.cpp
  xyzplugin.cpp
)

add_chemkit_plugin(xyz ${SOURCES})
target_link_libraries(xyz chemkit ${QT_LIBRARIES})
