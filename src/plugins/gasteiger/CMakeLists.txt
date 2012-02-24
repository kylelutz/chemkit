find_package(Chemkit REQUIRED)
include_directories(${CHEMKIT_INCLUDE_DIRS})

set(SOURCES
  gasteigerpartialchargepredictor.cpp
  gasteigerplugin.cpp
)

add_chemkit_plugin(gasteiger ${SOURCES})
target_link_libraries(gasteiger ${CHEMKIT_LIBRARIES})
