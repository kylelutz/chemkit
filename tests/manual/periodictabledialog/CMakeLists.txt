find_package(Qt4 4.6 COMPONENTS QtCore QtGui QtTest REQUIRED)
set(QT_USE_QTTEST TRUE)
include(${QT_USE_FILE})

include_directories(../../../include)

qt4_wrap_cpp(MOC_SOURCES periodictabledialogtest.h)
add_executable(periodictabledialogtest periodictabledialogtest.cpp ${MOC_SOURCES})
target_link_libraries(periodictabledialogtest chemkit chemkit-widgets ${QT_LIBRARIES})
