qt4_wrap_cpp(MOC_SOURCES rof5test.h)
add_executable(rof5test rof5test.cpp ${MOC_SOURCES})
target_link_libraries(rof5test chemkit ${QT_LIBRARIES})
add_chemkit_test(plugins.Rof5 rof5test)
