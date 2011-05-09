find_package(Qt4 4.6 COMPONENTS QtCore REQUIRED)
set(QT_DONT_USE_QTGUI TRUE)
include(${QT_USE_FILE})

find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

set(SOURCES
  xyzfileformat.cpp
  xyzplugin.cpp
)

set(MOC_HEADERS
  xyzplugin.h
)

qt4_wrap_cpp(MOC_SOURCES ${MOC_HEADERS})

add_chemkit_plugin(xyz ${SOURCES} ${MOC_SOURCES})
target_link_libraries(xyz chemkit ${QT_LIBRARIES})
