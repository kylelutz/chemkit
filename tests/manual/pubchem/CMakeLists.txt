if(NOT ${CHEMKIT_WITH_WEB})
  return()
endif()

find_package(Qt4 4.6 COMPONENTS QtCore QtNetwork QtTest REQUIRED)
set(QT_USE_QTNETWORK TRUE)
set(QT_USE_QTTEST TRUE)
include(${QT_USE_FILE})

find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

include_directories(../../../include)

qt4_wrap_cpp(MOC_SOURCES pubchemtest.h)
add_executable(pubchemtest pubchemtest.cpp ${MOC_SOURCES})
target_link_libraries(pubchemtest chemkit chemkit-io chemkit-web ${QT_LIBRARIES})
