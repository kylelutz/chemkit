find_package(Qt4 4.6 COMPONENTS QtCore REQUIRED)
set(QT_DONT_USE_QTGUI TRUE)
include(${QT_USE_FILE})

find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

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

set(MOC_HEADERS
  mmffplugin.h
)

qt4_wrap_cpp(MOC_SOURCES ${MOC_HEADERS})

add_chemkit_plugin(mmff ${SOURCES} ${MOC_SOURCES})
target_link_libraries(mmff chemkit ${QT_LIBRARIES})
