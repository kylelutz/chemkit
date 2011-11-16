qt4_wrap_cpp(MOC_SOURCES conformertest.h)
add_executable(conformertest conformertest.cpp ${MOC_SOURCES})
target_link_libraries(conformertest chemkit ${QT_LIBRARIES})
add_chemkit_test(chemkit.Conformer conformertest)
