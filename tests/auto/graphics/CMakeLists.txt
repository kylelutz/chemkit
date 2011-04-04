find_package(Qt4 4.6 COMPONENTS QtCore QtGui QtOpenGL QtTest REQUIRED)
set(QT_USE_QTOPENGL TRUE)
set(QT_USE_QTTEST TRUE)
include(${QT_USE_FILE})

include_directories(../../../include)

add_subdirectory(graphicscamera)
add_subdirectory(graphicsray)
add_subdirectory(graphicsscene)
add_subdirectory(graphicstransform)
