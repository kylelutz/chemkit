if(NOT ${CHEMKIT_WITH_IO})
  return()
endif()

qt4_wrap_cpp(MOC_SOURCES fhztest.h)
add_executable(fhztest fhztest.cpp ${MOC_SOURCES})
target_link_libraries(fhztest chemkit chemkit-io ${QT_LIBRARIES})
add_chemkit_test(fhz fhztest)
