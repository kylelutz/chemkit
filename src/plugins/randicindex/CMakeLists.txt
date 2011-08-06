find_package(Qt4 4.6 COMPONENTS QtCore REQUIRED)
set(QT_DONT_USE_QTGUI TRUE)
include(${QT_USE_FILE})

find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

set(SOURCES
  randicindexdescriptor.cpp
  randicindexplugin.cpp
)

add_chemkit_plugin(randicindex ${SOURCES})
target_link_libraries(randicindex chemkit ${QT_LIBRARIES})
