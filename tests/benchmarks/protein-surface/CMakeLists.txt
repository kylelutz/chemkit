if(NOT ${CHEMKIT_WITH_IO})
  return()
endif()

find_package(Chemkit COMPONENTS io)
include_directories(${CHEMKIT_INCLUDE_DIRS})

find_package(Qt4 4.6 COMPONENTS QtCore QtTest REQUIRED)
set(QT_DONT_USE_QTGUI TRUE)
set(QT_USE_QTTEST TRUE)
include(${QT_USE_FILE})

qt4_wrap_cpp(MOC_SOURCES proteinsurfacebenchmark.h)
add_executable(proteinsurfacebenchmark proteinsurfacebenchmark.cpp ${MOC_SOURCES})
target_link_libraries(proteinsurfacebenchmark ${CHEMKIT_LIBRARIES} ${QT_LIBRARIES})
