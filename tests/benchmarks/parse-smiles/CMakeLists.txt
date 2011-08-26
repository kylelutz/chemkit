if(NOT ${CHEMKIT_WITH_IO})
  return()
endif()

find_package(Qt4 4.6 COMPONENTS QtCore QtTest REQUIRED)
set(QT_DONT_USE_QTGUI TRUE)
set(QT_USE_QTTEST TRUE)
include(${QT_USE_FILE})

find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

include_directories(../../../include)

qt4_wrap_cpp(MOC_SOURCES parsesmilesbenchmark.h)
add_executable(parsesmilesbenchmark parsesmilesbenchmark.cpp ${MOC_SOURCES})
target_link_libraries(parsesmilesbenchmark chemkit chemkit-io ${QT_LIBRARIES})

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/pubchem-smiles.smi DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
