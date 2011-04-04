find_package(Qt4 4.6 COMPONENTS QtCore QtTest REQUIRED)
set(QT_DONT_USE_QTGUI TRUE)
set(QT_USE_QTTEST TRUE)
include(${QT_USE_FILE})

include_directories(../../../include)

qt4_wrap_cpp(MOC_SOURCES mmffenergybenchmark.h)
add_executable(mmffenergybenchmark mmffenergybenchmark.cpp ${MOC_SOURCES})
target_link_libraries(mmffenergybenchmark chemkit ${QT_LIBRARIES})
