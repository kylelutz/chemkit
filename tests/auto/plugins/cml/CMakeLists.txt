qt4_wrap_cpp(MOC_SOURCES cmltest.h)
add_executable(cmltest cmltest.cpp ${MOC_SOURCES})
target_link_libraries(cmltest chemkit ${QT_LIBRARIES})
add_chemkit_test(cml cmltest)
