qt4_wrap_cpp(MOC_SOURCES isotopetest.h)
add_executable(isotopetest isotopetest.cpp ${MOC_SOURCES})
target_link_libraries(isotopetest chemkit ${QT_LIBRARIES})
add_chemkit_test(chemkit.Isotope isotopetest)
