qt4_wrap_cpp(MOC_SOURCES forcefieldtest.h)
add_executable(forcefieldtest forcefieldtest.cpp mockforcefield.cpp ${MOC_SOURCES})
target_link_libraries(forcefieldtest chemkit chemkit-md ${QT_LIBRARIES})
add_chemkit_test(md.ForceField forcefieldtest)
