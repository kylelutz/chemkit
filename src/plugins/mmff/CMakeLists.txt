find_package(Chemkit COMPONENTS md REQUIRED)
include_directories(${CHEMKIT_INCLUDE_DIRS})

find_package(Qt4 4.6 COMPONENTS QtCore REQUIRED)
set(QT_DONT_USE_QTGUI TRUE)
include(${QT_USE_FILE})

set(SOURCES
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
target_link_libraries(mmff ${CHEMKIT_LIBRARIES} ${QT_LIBRARIES})
