if(NOT ${CHEMKIT_WITH_IO})
  return()
endif()

qt4_wrap_cpp(MOC_SOURCES pqrtest.h)
add_executable(pqrtest pqrtest.cpp ${MOC_SOURCES})
target_link_libraries(pqrtest chemkit chemkit-io ${QT_LIBRARIES})
add_chemkit_test(plugins.Pqr pqrtest)
