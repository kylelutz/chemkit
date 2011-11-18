find_package(Chemkit)
include_directories(${CHEMKIT_INCLUDE_DIRS})

qt4_wrap_cpp(MOC_SOURCES translatetest.h)
add_executable(translatetest translatetest.cpp ${MOC_SOURCES})
target_link_libraries(translatetest ${CHEMKIT_LIBRARIES} ${QT_LIBRARIES})
add_chemkit_test(apps.Translate translatetest)
