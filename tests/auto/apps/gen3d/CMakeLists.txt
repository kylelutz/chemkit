find_package(Chemkit COMPONENTS io)
include_directories(${CHEMKIT_INCLUDE_DIRS})

qt4_wrap_cpp(MOC_SOURCES gen3dtest.h)
add_executable(gen3dtest gen3dtest.cpp ${MOC_SOURCES})
target_link_libraries(gen3dtest ${CHEMKIT_LIBRARIES} ${QT_LIBRARIES})
add_chemkit_test(apps.Gen3d gen3dtest)
