qt4_wrap_cpp(MOC_SOURCES translatetest.h)
add_executable(translatetest translatetest.cpp ${MOC_SOURCES})
target_link_libraries(translatetest chemkit ${QT_LIBRARIES})
add_chemkit_test(translate translatetest)
