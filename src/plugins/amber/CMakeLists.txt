find_package(Qt4 4.6 COMPONENTS QtCore REQUIRED)
set(QT_DONT_USE_QTGUI TRUE)
include(${QT_USE_FILE})

set(SOURCES
  ambercalculation.cpp
  amberforcefield.cpp
  amberparameters.cpp
  amberplugin.cpp
)

set(MOC_HEADERS
  amberplugin.h
)

qt4_wrap_cpp(MOC_SOURCES ${MOC_HEADERS})

add_chemkit_plugin(amber ${SOURCES} ${MOC_SOURCES})
target_link_libraries(amber chemkit ${QT_LIBRARIES})
