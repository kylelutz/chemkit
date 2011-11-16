qt4_wrap_cpp(MOC_SOURCES atomtest.h)
add_executable(atomtest atomtest.cpp ${MOC_SOURCES})
target_link_libraries(atomtest chemkit ${QT_LIBRARIES})
add_chemkit_test(chemkit.Atom atomtest)
