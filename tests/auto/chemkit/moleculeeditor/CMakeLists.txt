qt4_wrap_cpp(MOC_SOURCES moleculeeditortest.h)
add_executable(moleculeeditortest moleculeeditortest.cpp ${MOC_SOURCES})
target_link_libraries(moleculeeditortest chemkit ${QT_LIBRARIES})
add_chemkit_test(chemkit.MoleculeEditor moleculeeditortest)
