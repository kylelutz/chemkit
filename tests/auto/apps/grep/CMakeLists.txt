find_package(Chemkit COMPONENTS io)
include_directories(${CHEMKIT_INCLUDE_DIRS})

qt4_wrap_cpp(MOC_SOURCES greptest.h)
add_executable(greptest greptest.cpp ${MOC_SOURCES})
target_link_libraries(greptest ${CHEMKIT_LIBRARIES} ${QT_LIBRARIES})
add_chemkit_test(apps.Grep greptest)
