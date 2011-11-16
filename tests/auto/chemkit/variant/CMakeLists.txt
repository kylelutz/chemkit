qt4_wrap_cpp(MOC_SOURCES varianttest.h)
add_executable(varianttest varianttest.cpp ${MOC_SOURCES})
target_link_libraries(varianttest chemkit ${QT_LIBRARIES})
add_chemkit_test(chemkit.Variant varianttest)
