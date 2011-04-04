qt4_wrap_cpp(MOC_SOURCES staticmatrixtest.h)
add_executable(staticmatrixtest staticmatrixtest.cpp ${MOC_SOURCES})
target_link_libraries(staticmatrixtest chemkit ${QT_LIBRARIES})
add_chemkit_test(staticmatrix staticmatrixtest)
