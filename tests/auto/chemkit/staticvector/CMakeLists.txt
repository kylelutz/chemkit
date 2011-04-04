qt4_wrap_cpp(MOC_SOURCES staticvectortest.h)
add_executable(staticvectortest staticvectortest.cpp ${MOC_SOURCES})
target_link_libraries(staticvectortest chemkit ${QT_LIBRARIES})
add_chemkit_test(staticvector staticvectortest)
