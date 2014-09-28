find_package(Chemkit COMPONENTS graphics)
include_directories(${CHEMKIT_INCLUDE_DIRS})

find_package(Qt5 COMPONENTS Widgets OpenGL Designer REQUIRED)
set(QT_USE_QTOPENGL TRUE)
set(QT_USE_QTDESIGNER TRUE)
include_directories(${Qt5Widgets_INCLUDE_DIRS})
include_directories(${Qt5OpenGL_INCLUDE_DIRS})
include_directories(${Qt5Designer_INCLUDE_DIRS})

set(SOURCES graphicsviewdesignerplugin.cpp)
set(MOC_HEADERS graphicsviewdesignerplugin.h)

qt5_wrap_cpp(MOC_SOURCES ${MOC_HEADERS})

add_library(graphicsviewdesignerplugin SHARED ${SOURCES} ${MOC_SOURCES})
target_link_libraries(graphicsviewdesignerplugin Qt5::Widgets Qt5::OpenGL Qt5::Designer ${CHEMKIT_LIBRARIES})
set_target_properties(graphicsviewdesignerplugin PROPERTIES OUTPUT_NAME chemkitgraphicsviewplugin)
#install(TARGET graphicsviewdesignerplugin DESTINATION ${QT_INSTALL_PLUGINS}/designer)
