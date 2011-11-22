qt4_wrap_cpp(MOC_SOURCES graphdescriptorstest.h)
add_executable(graphdescriptorstest graphdescriptorstest.cpp ${MOC_SOURCES})
target_link_libraries(graphdescriptorstest chemkit ${QT_LIBRARIES})
add_chemkit_test(plugins.GraphDescriptors graphdescriptorstest)
