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

set(MOC_HEADERS
  oplsplugin.h
)

qt4_wrap_cpp(MOC_SOURCES ${MOC_HEADERS})

add_chemkit_plugin(opls ${SOURCES} ${MOC_SOURCES})
target_link_libraries(opls chemkit ${QT_LIBRARIES})
