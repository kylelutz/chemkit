if(NOT ${CHEMKIT_WITH_IO})
  return()
endif()

find_package(Chemkit COMPONENTS io REQUIRED)
include_directories(${CHEMKIT_INCLUDE_DIRS})

set(SOURCES
  chemjsonfileformat.cpp
  chemjsonplugin.cpp
)

aux_source_directory(../../3rdparty/jsoncpp/ JSONCPP_SOURCES)

add_chemkit_plugin(chemjson ${SOURCES} ${JSONCPP_SOURCES})
target_link_libraries(chemjson ${CHEMKIT_LIBRARIES})
