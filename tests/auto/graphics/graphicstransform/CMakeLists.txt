qt4_wrap_cpp(MOC_SOURCES graphicstransformtest.h)
add_executable(graphicstransformtest graphicstransformtest.cpp ${MOC_SOURCES})
target_link_libraries(graphicstransformtest chemkit chemkit-graphics ${QT_LIBRARIES})
add_chemkit_test(graphicstransform graphicstransformtest)
