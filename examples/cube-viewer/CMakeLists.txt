find_package(Qt4 4.6 COMPONENTS QtCore QtGui QtOpenGL REQUIRED)
set(QT_USE_QTOPENGL TRUE)
include(${QT_USE_FILE})

find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

include_directories(${CMAKE_CURRENT_BINARY_DIR})
include_directories(../../include)

set(SOURCES cubeviewerexample.cpp)
set(MOC_HEADERS cubeviewerexample.h)
set(UI_FORMS cubeviewerexample.ui)

qt4_wrap_ui(UI_SOURCES ${UI_FORMS})
qt4_wrap_cpp(MOC_SOURCES ${MOC_HEADERS})

add_executable(cube-viewer ${SOURCES} ${MOC_SOURCES} ${UI_SOURCES})
target_link_libraries(cube-viewer chemkit chemkit-io chemkit-graphics ${QT_LIBRARIES})
