qt4_wrap_cpp(MOC_SOURCES surfacedescriptorstest.h)
add_executable(surfacedescriptorstest surfacedescriptorstest.cpp ${MOC_SOURCES})
target_link_libraries(surfacedescriptorstest chemkit chemkit-io ${QT_LIBRARIES})
add_chemkit_test(plugins.SurfaceDescriptors surfacedescriptorstest)
