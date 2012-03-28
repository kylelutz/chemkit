qt4_wrap_cpp(MOC_SOURCES elementtyperstest.h)
add_executable(elementtyperstest elementtyperstest.cpp ${MOC_SOURCES})
target_link_libraries(elementtyperstest chemkit ${QT_LIBRARIES})
add_chemkit_test(plugins.ElementTypers elementtyperstest)
