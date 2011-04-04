qt4_wrap_cpp(MOC_SOURCES coordinatestest.h)
add_executable(coordinatestest coordinatestest.cpp ${MOC_SOURCES})
target_link_libraries(coordinatestest chemkit ${QT_LIBRARIES})
add_chemkit_test(coordinates coordinatestest)
