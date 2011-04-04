find_package(Qt4 4.6 COMPONENTS QtCore REQUIRED)
set(QT_DONT_USE_QTGUI TRUE)
include(${QT_USE_FILE})

include_directories(../../../include)

add_chemkit_executable(convert convert.cpp)
target_link_libraries(convert chemkit ${QT_LIBRARIES})
