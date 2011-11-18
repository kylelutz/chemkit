if(NOT ${CHEMKIT_WITH_IO} OR NOT ${CHEMKIT_WITH_MD})
  return()
endif()

find_package(Chemkit COMPONENTS io md)
include_directories(${CHEMKIT_INCLUDE_DIRS})

find_package(Qt4 4.6 COMPONENTS QtCore REQUIRED)
set(QT_DONT_USE_QTGUI TRUE)
include(${QT_USE_FILE})

include_directories(${CMAKE_CURRENT_BINARY_DIR})

set(SOURCES uffenergyexample.cpp)

add_executable(uff-energy ${SOURCES})
target_link_libraries(uff-energy ${CHEMKIT_LIBRARIES} ${QT_LIBRARIES})
