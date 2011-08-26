if(NOT ${CHEMKIT_WITH_IO})
  return()
endif()

qt4_wrap_cpp(MOC_SOURCES xyztest.h)
add_executable(xyztest xyztest.cpp ${MOC_SOURCES})
target_link_libraries(xyztest chemkit chemkit-io ${QT_LIBRARIES})
add_chemkit_test(xyz xyztest)
