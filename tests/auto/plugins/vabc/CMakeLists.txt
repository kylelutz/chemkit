qt4_wrap_cpp(MOC_SOURCES vabctest.h)
add_executable(vabctest vabctest.cpp ${MOC_SOURCES})
target_link_libraries(vabctest chemkit ${QT_LIBRARIES})
add_chemkit_test(plugins.Vabc vabctest)
