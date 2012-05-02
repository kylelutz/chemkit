if(NOT ${CHEMKIT_WITH_IO})
  return()
endif()

qt4_wrap_cpp(MOC_SOURCES chemjsontest.h)
add_executable(chemjsontest chemjsontest.cpp ${MOC_SOURCES})
target_link_libraries(chemjsontest chemkit chemkit-io ${QT_LIBRARIES})
add_chemkit_test(plugins.ChemJson chemjsontest)
