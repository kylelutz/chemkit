find_package(Chemkit COMPONENTS io OPTIONAL)
include_directories(${CHEMKIT_INCLUDE_DIRS})

find_package(Qt4 4.6 COMPONENTS QtCore REQUIRED)
set(QT_DONT_USE_QTGUI TRUE)
include(${QT_USE_FILE})

set(SOURCES
  sybylatomtyper.cpp
  sybylplugin.cpp
)

if(${CHEMKIT_WITH_IO})
  set(SOURCES ${SOURCES} mol2fileformat.cpp)
endif()

add_chemkit_plugin(sybyl ${SOURCES})
target_link_libraries(sybyl ${CHEMKIT_LIBRARIES} ${QT_LIBRARIES})
