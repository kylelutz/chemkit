qt4_wrap_cpp(MOC_SOURCES greptest.h)
add_executable(greptest greptest.cpp ${MOC_SOURCES})
target_link_libraries(greptest chemkit ${QT_LIBRARIES})
add_chemkit_test(grep greptest)
