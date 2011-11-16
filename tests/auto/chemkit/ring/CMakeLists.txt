qt4_wrap_cpp(MOC_SOURCES ringtest.h)
add_executable(ringtest ringtest.cpp ${MOC_SOURCES})
target_link_libraries(ringtest chemkit ${QT_LIBRARIES})
add_chemkit_test(chemkit.Ring ringtest)
