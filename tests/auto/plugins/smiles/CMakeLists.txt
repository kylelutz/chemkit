qt4_wrap_cpp(MOC_SOURCES smilestest.h)
add_executable(smilestest smilestest.cpp ${MOC_SOURCES})
target_link_libraries(smilestest chemkit ${QT_LIBRARIES})
add_chemkit_test(smiles smilestest)
