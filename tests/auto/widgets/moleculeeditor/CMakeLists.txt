qt4_wrap_cpp(MOC_SOURCES moleculeeditortest.h)
add_executable(moleculeeditortest moleculeeditortest.cpp ${MOC_SOURCES})
target_link_libraries(moleculeeditortest chemkit chemkit-widgets ${QT_LIBRARIES})
add_chemkit_test(widgets.MoleculeEditor moleculeeditortest)
