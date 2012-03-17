qt4_wrap_cpp(MOC_SOURCES mlogptest.h)
add_executable(mlogptest mlogptest.cpp ${MOC_SOURCES})
target_link_libraries(mlogptest chemkit ${QT_LIBRARIES})
add_chemkit_test(plugins.MlogP mlogptest)
