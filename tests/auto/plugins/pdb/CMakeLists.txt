qt4_wrap_cpp(MOC_SOURCES pdbtest.h)
add_executable(pdbtest pdbtest.cpp ${MOC_SOURCES})
target_link_libraries(pdbtest chemkit chemkit-io ${QT_LIBRARIES})
add_chemkit_test(pdb pdbtest)
