find_package(Chemkit COMPONENTS md REQUIRED)
include_directories(${CHEMKIT_INCLUDE_DIRS})

find_package(Qt4 4.6 COMPONENTS QtCore REQUIRED)
set(QT_DONT_USE_QTGUI TRUE)
include(${QT_USE_FILE})

set(SOURCES
  oplsatomtyper.cpp
  oplscalculation.cpp
  oplsforcefield.cpp
  oplsparameters.cpp
  oplsplugin.cpp
)

add_chemkit_plugin(opls ${SOURCES})
target_link_libraries(opls ${CHEMKIT_LIBRARIES} ${QT_LIBRARIES})
