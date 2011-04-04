qt4_wrap_cpp(MOC_SOURCES fragmenttest.h)
add_executable(fragmenttest fragmenttest.cpp ${MOC_SOURCES})
target_link_libraries(fragmenttest chemkit ${QT_LIBRARIES})
add_chemkit_test(fragment fragmenttest)
