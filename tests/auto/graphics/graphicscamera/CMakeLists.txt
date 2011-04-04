qt4_wrap_cpp(MOC_SOURCES graphicscameratest.h)
add_executable(graphicscameratest graphicscameratest.cpp ${MOC_SOURCES})
target_link_libraries(graphicscameratest chemkit chemkit-graphics ${QT_LIBRARIES})
add_chemkit_test(graphicscamera graphicscameratest)
