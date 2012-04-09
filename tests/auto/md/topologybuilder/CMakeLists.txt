qt4_wrap_cpp(MOC_SOURCES topologybuildertest.h)
add_executable(topologybuildertest topologybuildertest.cpp ${MOC_SOURCES})
target_link_libraries(topologybuildertest chemkit chemkit-md ${QT_LIBRARIES})
add_chemkit_test(md.TopologyBuilder topologybuildertest)
