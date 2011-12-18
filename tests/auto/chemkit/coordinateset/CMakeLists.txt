qt4_wrap_cpp(MOC_SOURCES coordinatesettest.h)
add_executable(coordinatesettest coordinatesettest.cpp ${MOC_SOURCES})
target_link_libraries(coordinatesettest chemkit ${QT_LIBRARIES})
add_chemkit_test(chemkit.CoordinateSet coordinatesettest)
