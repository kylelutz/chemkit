find_package(Chemkit COMPONENTS md REQUIRED)
include_directories(${CHEMKIT_INCLUDE_DIRS})

set(SOURCES
  mmffaromaticitymodel.cpp
  mmffatom.cpp
  mmffatomtyper.cpp
  mmffcalculation.cpp
  mmffforcefield.cpp
  mmffparameters.cpp
  mmffparametersdata.cpp
  mmffpartialchargepredictor.cpp
  mmffplugin.cpp
)

add_chemkit_plugin(mmff ${SOURCES})
target_link_libraries(mmff ${CHEMKIT_LIBRARIES})
