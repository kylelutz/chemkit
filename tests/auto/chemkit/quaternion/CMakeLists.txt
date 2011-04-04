qt4_wrap_cpp(MOC_SOURCES quaterniontest.h)
add_executable(quaterniontest quaterniontest.cpp ${MOC_SOURCES})
target_link_libraries(quaterniontest chemkit ${QT_LIBRARIES})
add_chemkit_test(quaternion quaterniontest)
