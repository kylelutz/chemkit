qt4_wrap_cpp(MOC_SOURCES converttest.h)
add_executable(converttest converttest.cpp ${MOC_SOURCES})
target_link_libraries(converttest chemkit ${QT_LIBRARIES})
add_chemkit_test(convert converttest)
