if(NOT ${CHEMKIT_WITH_IO})
  return()
endif()

find_package(Chemkit COMPONENTS io)
include_directories(${CHEMKIT_INCLUDE_DIRS})

find_package(Qt4 4.6 COMPONENTS QtCore QtTest REQUIRED)
set(QT_DONT_USE_QTGUI TRUE)
set(QT_USE_QTTEST TRUE)
include(${QT_USE_FILE})

qt4_wrap_cpp(MOC_SOURCES molecularmassesbenchmark.h)
add_executable(molecularmassesbenchmark molecularmassesbenchmark.cpp ${MOC_SOURCES})
target_link_libraries(molecularmassesbenchmark ${CHEMKIT_LIBRARIES} ${QT_LIBRARIES})

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/pubchem_sample_33.sdf DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
