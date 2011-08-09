qt4_wrap_cpp(MOC_SOURCES cubetest.h)
add_executable(cubetest cubetest.cpp ${MOC_SOURCES})
target_link_libraries(cubetest chemkit chemkit-io ${QT_LIBRARIES})
add_chemkit_test(cube cubetest)
