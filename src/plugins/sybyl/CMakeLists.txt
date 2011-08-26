find_package(Qt4 4.6 COMPONENTS QtCore REQUIRED)
set(QT_DONT_USE_QTGUI TRUE)
include(${QT_USE_FILE})

find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

set(SOURCES
  sybylatomtyper.cpp
  sybylplugin.cpp
)

if(${CHEMKIT_WITH_IO})
  set(SOURCES ${SOURCES} mol2fileformat.cpp)
endif()

add_chemkit_plugin(sybyl ${SOURCES})
target_link_libraries(sybyl chemkit ${QT_LIBRARIES})

if(${CHEMKIT_WITH_IO})
  target_link_libraries(sybyl chemkit-io)
endif()
