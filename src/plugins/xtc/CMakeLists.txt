if(NOT ${CHEMKIT_WITH_MD_IO})
  return()
endif()

find_package(Boost COMPONENTS system filesystem REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

find_package(Chemkit COMPONENTS io md md-io REQUIRED)
include_directories(${CHEMKIT_INCLUDE_DIRS})

set(SOURCES
  xtcfileformat.cpp
  xtcplugin.cpp
)

aux_source_directory(../../3rdparty/xdrf/ XDRF_SOURCES)

add_chemkit_plugin(xtc ${SOURCES} ${XDRF_SOURCES})
target_link_libraries(xtc ${CHEMKIT_LIBRARIES} ${Boost_LIBRARIES})
