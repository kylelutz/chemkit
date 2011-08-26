if(NOT ${CHEMKIT_WITH_IO})
  return()
endif()

qt4_wrap_cpp(MOC_SOURCES mdltest.h)
add_executable(mdltest mdltest.cpp ${MOC_SOURCES})
target_link_libraries(mdltest chemkit chemkit-io ${QT_LIBRARIES})
add_chemkit_test(mdl mdltest)
