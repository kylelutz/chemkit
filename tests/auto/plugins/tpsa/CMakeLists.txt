qt4_wrap_cpp(MOC_SOURCES tpsatest.h)
add_executable(tpsatest tpsatest.cpp ${MOC_SOURCES})
target_link_libraries(tpsatest chemkit ${QT_LIBRARIES})
add_chemkit_test(plugins.Tpsa tpsatest)
