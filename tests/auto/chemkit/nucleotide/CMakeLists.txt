qt4_wrap_cpp(MOC_SOURCES nucleotidetest.h)
add_executable(nucleotidetest nucleotidetest.cpp ${MOC_SOURCES})
target_link_libraries(nucleotidetest chemkit ${QT_LIBRARIES})
add_chemkit_test(chemkit.Nucleotide nucleotidetest)
