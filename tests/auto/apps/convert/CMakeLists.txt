find_package(Chemkit COMPONENTS io)
include_directories(${CHEMKIT_INCLUDE_DIRS})

qt4_wrap_cpp(MOC_SOURCES converttest.h)
add_executable(converttest converttest.cpp ${MOC_SOURCES})
target_link_libraries(converttest ${CHEMKIT_LIBRARIES} ${QT_LIBRARIES})
add_chemkit_test(apps.Convert converttest)
