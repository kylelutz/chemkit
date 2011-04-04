qt4_wrap_cpp(MOC_SOURCES moleculetest.h)
add_executable(moleculetest moleculetest.cpp ${MOC_SOURCES})
target_link_libraries(moleculetest chemkit ${QT_LIBRARIES})
add_chemkit_test(molecule moleculetest)
