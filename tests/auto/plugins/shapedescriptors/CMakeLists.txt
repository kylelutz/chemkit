qt4_wrap_cpp(MOC_SOURCES shapedescriptorstest.h)
add_executable(shapedescriptorstest shapedescriptorstest.cpp ${MOC_SOURCES})
target_link_libraries(shapedescriptorstest chemkit chemkit-io ${QT_LIBRARIES})
add_chemkit_test(plugins.ShapeDescriptors shapedescriptorstest)
