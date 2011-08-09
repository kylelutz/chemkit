qt4_wrap_cpp(MOC_SOURCES molecularsurfacetest.h)
add_executable(molecularsurfacetest molecularsurfacetest.cpp ${MOC_SOURCES})
target_link_libraries(molecularsurfacetest chemkit chemkit-io ${QT_LIBRARIES})
add_chemkit_test(molecularsurface molecularsurfacetest)
