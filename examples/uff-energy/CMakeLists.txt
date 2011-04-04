find_package(Qt4 4.6 COMPONENTS QtCore REQUIRED)
set(QT_DONT_USE_QTGUI TRUE)
include(${QT_USE_FILE})

include_directories(${CMAKE_CURRENT_BINARY_DIR})
include_directories(../../include)

set(SOURCES uffenergyexample.cpp)

add_executable(uff-energy ${SOURCES})
target_link_libraries(uff-energy chemkit ${QT_LIBRARIES})
