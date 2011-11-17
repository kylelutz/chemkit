find_package(Chemkit COMPONENTS md OPTIONAL)
include_directories(${CHEMKIT_INCLUDE_DIRS})

set(SOURCES
  uffatomtyper.cpp
  uffplugin.cpp
)

if(${CHEMKIT_WITH_MD})
  set(SOURCES ${SOURCES}
    uffcalculation.cpp
    uffforcefield.cpp
    uffparameters.cpp
  )
endif()

add_chemkit_plugin(uff ${SOURCES})
target_link_libraries(uff ${CHEMKIT_LIBRARIES})
