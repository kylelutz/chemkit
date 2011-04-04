qt4_wrap_cpp(MOC_SOURCES point3test.h)
add_executable(point3test point3test.cpp ${MOC_SOURCES})
target_link_libraries(point3test chemkit ${QT_LIBRARIES})
add_chemkit_test(point3 point3test)
