if(NOT ${CHEMKIT_WITH_IO})
  return()
endif()

qt4_wrap_cpp(MOC_SOURCES fpstest.h)
add_executable(fpstest fpstest.cpp ${MOC_SOURCES})
target_link_libraries(fpstest chemkit chemkit-io ${QT_LIBRARIES})
add_chemkit_test(plugins.Fps fpstest)
