qt4_wrap_cpp(MOC_SOURCES polymertest.h)
add_executable(polymertest polymertest.cpp ${MOC_SOURCES})
target_link_libraries(polymertest chemkit ${QT_LIBRARIES})
add_chemkit_test(chemkit.Polymer polymertest)
