find_package(Chemkit REQUIRED)
include_directories(${CHEMKIT_INCLUDE_DIRS})

find_package(Qt4 4.6 COMPONENTS QtCore REQUIRED)
set(QT_DONT_USE_QTGUI TRUE)
include(${QT_USE_FILE})

set(SOURCES
  gasteigerpartialchargepredictor.cpp
  gasteigerplugin.cpp
)

add_chemkit_plugin(gasteiger ${SOURCES})
target_link_libraries(gasteiger ${CHEMKIT_LIBRARIES} ${QT_LIBRARIES})
