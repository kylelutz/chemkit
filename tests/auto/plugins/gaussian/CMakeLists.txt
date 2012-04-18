if(NOT ${CHEMKIT_WITH_IO})
  return()
endif()

qt4_wrap_cpp(MOC_SOURCES gaussiantest.h)
add_executable(gaussiantest gaussiantest.cpp ${MOC_SOURCES})
target_link_libraries(gaussiantest chemkit chemkit-io ${QT_LIBRARIES})
add_chemkit_test(plugins.Gaussian gaussiantest)
