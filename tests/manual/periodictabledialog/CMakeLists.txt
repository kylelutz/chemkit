if(NOT ${CHEMKIT_WITH_GUI})
  return()
endif()

find_package(Chemkit COMPONENTS gui)
include_directories(${CHEMKIT_INCLUDE_DIRS})

find_package(Qt4 4.6 COMPONENTS QtCore QtGui QtTest REQUIRED)
set(QT_USE_QTTEST TRUE)
include(${QT_USE_FILE})

qt4_wrap_cpp(MOC_SOURCES periodictabledialogtest.h)
add_executable(periodictabledialogtest periodictabledialogtest.cpp ${MOC_SOURCES})
target_link_libraries(periodictabledialogtest ${CHEMKIT_LIBRARIES} ${QT_LIBRARIES})
