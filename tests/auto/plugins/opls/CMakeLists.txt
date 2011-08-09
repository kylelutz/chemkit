qt4_wrap_cpp(MOC_SOURCES oplstest.h)
add_executable(oplstest oplstest.cpp ${MOC_SOURCES})
target_link_libraries(oplstest chemkit chemkit-io ${QT_LIBRARIES})
add_chemkit_test(opls oplstest)
