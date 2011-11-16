if(NOT ${CHEMKIT_WITH_IO})
  return()
endif()

qt4_wrap_cpp(MOC_SOURCES txyztest.h)
add_executable(txyztest txyztest.cpp ${MOC_SOURCES})
target_link_libraries(txyztest chemkit chemkit-io ${QT_LIBRARIES})
add_chemkit_test(plugins.Txyz txyztest)
