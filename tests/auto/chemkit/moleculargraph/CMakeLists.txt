qt4_wrap_cpp(MOC_SOURCES moleculargraphtest.h)
add_executable(moleculargraphtest moleculargraphtest.cpp ${MOC_SOURCES})
target_link_libraries(moleculargraphtest chemkit ${QT_LIBRARIES})
add_chemkit_test(chemkit.MolecularGraph moleculargraphtest)
