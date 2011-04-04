qt4_wrap_cpp(MOC_SOURCES vector3test.h)
add_executable(vector3test vector3test.cpp ${MOC_SOURCES})
target_link_libraries(vector3test chemkit ${QT_LIBRARIES})
add_chemkit_test(vector3 vector3test)
