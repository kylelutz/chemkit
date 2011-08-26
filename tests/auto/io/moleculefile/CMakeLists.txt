qt4_wrap_cpp(MOC_SOURCES moleculefiletest.h)
add_executable(moleculefiletest moleculefiletest.cpp ${MOC_SOURCES})
target_link_libraries(moleculefiletest chemkit chemkit-io ${QT_LIBRARIES})
add_chemkit_test(moleculefile moleculefiletest)
