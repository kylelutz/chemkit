find_package(Chemkit COMPONENTS graphics)
include_directories(${CHEMKIT_INCLUDE_DIRS})

find_package(Qt4 4.6 COMPONENTS QtCore QtGui QtOpenGL QtDesigner REQUIRED)
set(QT_USE_QTOPENGL TRUE)
set(QT_USE_QTDESIGNER TRUE)
include(${QT_USE_FILE})

set(SOURCES graphicsviewdesignerplugin.cpp)
set(MOC_HEADERS graphicsviewdesignerplugin.h)

qt4_wrap_cpp(MOC_SOURCES ${MOC_HEADERS})

add_library(graphicsviewdesignerplugin SHARED ${SOURCES} ${MOC_SOURCES})
target_link_libraries(graphicsviewdesignerplugin ${CHEMKIT_LIBRARIES} ${QT_LIBRARIES})
set_target_properties(graphicsviewdesignerplugin PROPERTIES OUTPUT_NAME chemkitgraphicsviewplugin)
#install(TARGET graphicsviewdesignerplugin DESTINATION ${QT_INSTALL_PLUGINS}/designer)
