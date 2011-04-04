qt4_wrap_cpp(MOC_SOURCES elementtest.h)
add_executable(elementtest elementtest.cpp ${MOC_SOURCES})
target_link_libraries(elementtest chemkit ${QT_LIBRARIES})
add_chemkit_test(element elementtest)
