qt4_wrap_cpp(MOC_SOURCES mcdltest.h)
add_executable(mcdltest mcdltest.cpp ${MOC_SOURCES})
target_link_libraries(mcdltest chemkit ${QT_LIBRARIES})
add_chemkit_test(mcdl mcdltest)
