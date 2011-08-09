qt4_wrap_cpp(MOC_SOURCES pdbmltest.h)
add_executable(pdbmltest pdbmltest.cpp ${MOC_SOURCES})
target_link_libraries(pdbmltest chemkit chemkit-io ${QT_LIBRARIES})
add_chemkit_test(pdbml pdbmltest)
