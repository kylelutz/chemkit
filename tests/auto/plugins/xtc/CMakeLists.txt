if(NOT ${CHEMKIT_WITH_IO} OR NOT ${CHEMKIT_WITH_MD})
  return()
endif()

qt4_wrap_cpp(MOC_SOURCES xtctest.h)
add_executable(xtctest xtctest.cpp ${MOC_SOURCES})
target_link_libraries(xtctest chemkit chemkit-md ${QT_LIBRARIES})
add_chemkit_test(xtc xtctest)
