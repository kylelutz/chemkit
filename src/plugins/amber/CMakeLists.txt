find_package(Chemkit COMPONENTS md md-io REQUIRED)
include_directories(${CHEMKIT_INCLUDE_DIRS})

set(SOURCES
  amberatomtyper.cpp
  ambercalculation.cpp
  amberforcefield.cpp
  amberparameters.cpp
  amberplugin.cpp
)

if(${CHEMKIT_WITH_MD_IO})
  set(SOURCES ${SOURCES}
    mdcrdfileformat.cpp
  )
endif()

add_chemkit_plugin(amber ${SOURCES})
target_link_libraries(amber ${CHEMKIT_LIBRARIES})
