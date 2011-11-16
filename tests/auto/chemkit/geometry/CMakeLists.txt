qt4_wrap_cpp(MOC_SOURCES geometrytest.h)
add_executable(geometrytest geometrytest.cpp ${MOC_SOURCES})
target_link_libraries(geometrytest chemkit ${QT_LIBRARIES})
add_chemkit_test(chemkit.Geometry geometrytest)
