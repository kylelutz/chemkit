qt4_wrap_cpp(MOC_SOURCES topologytest.h)
add_executable(topologytest topologytest.cpp ${MOC_SOURCES})
target_link_libraries(topologytest chemkit chemkit-md ${QT_LIBRARIES})
add_chemkit_test(md.Topology topologytest)
