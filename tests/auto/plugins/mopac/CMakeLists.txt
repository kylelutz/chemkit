if(NOT ${CHEMKIT_WITH_IO})
  return()
endif()

qt4_wrap_cpp(MOC_SOURCES mopactest.h)
add_executable(mopactest mopactest.cpp ${MOC_SOURCES})
target_link_libraries(mopactest chemkit chemkit-io ${QT_LIBRARIES})
add_chemkit_test(plugins.Mopac mopactest)
