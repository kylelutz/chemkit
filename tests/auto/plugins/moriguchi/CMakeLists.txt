qt4_wrap_cpp(MOC_SOURCES moriguchitest.h)
add_executable(moriguchitest moriguchitest.cpp ${MOC_SOURCES})
target_link_libraries(moriguchitest chemkit ${QT_LIBRARIES})
add_chemkit_test(plugins.Moriguchi moriguchitest)
