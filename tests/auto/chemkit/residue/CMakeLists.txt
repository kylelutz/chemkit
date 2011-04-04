qt4_wrap_cpp(MOC_SOURCES residuetest.h)
add_executable(residuetest residuetest.cpp ${MOC_SOURCES})
target_link_libraries(residuetest chemkit ${QT_LIBRARIES})
add_chemkit_test(residue residuetest)
