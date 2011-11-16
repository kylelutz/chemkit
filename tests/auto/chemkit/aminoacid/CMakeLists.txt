qt4_wrap_cpp(MOC_SOURCES aminoacidtest.h)
add_executable(aminoacidtest aminoacidtest.cpp ${MOC_SOURCES})
target_link_libraries(aminoacidtest chemkit ${QT_LIBRARIES})
add_chemkit_test(chemkit.AminoAcid aminoacidtest)
