if(NOT ${CHEMKIT_WITH_MD_IO})
  return()
endif()

qt4_wrap_cpp(MOC_SOURCES gromacstest.h)
add_executable(gromacstest gromacstest.cpp ${MOC_SOURCES})
target_link_libraries(gromacstest chemkit chemkit-io chemkit-md chemkit-md-io ${QT_LIBRARIES})
add_chemkit_test(plugins.Gromacs gromacstest)
