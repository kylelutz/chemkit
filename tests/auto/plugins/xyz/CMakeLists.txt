qt4_wrap_cpp(MOC_SOURCES xyztest.h)
add_executable(xyztest xyztest.cpp ${MOC_SOURCES})
target_link_libraries(xyztest chemkit ${QT_LIBRARIES})
add_chemkit_test(xyz xyztest)
