qt4_wrap_cpp(MOC_SOURCES fp2test.h)
add_executable(fp2test fp2test.cpp ${MOC_SOURCES})
target_link_libraries(fp2test chemkit chemkit-io ${QT_LIBRARIES})
add_chemkit_test(plugins.Fp2 fp2test)
