qt4_wrap_cpp(MOC_SOURCES graphicsraytest.h)
add_executable(graphicsraytest graphicsraytest.cpp ${MOC_SOURCES})
target_link_libraries(graphicsraytest chemkit chemkit-graphics ${QT_LIBRARIES})
add_chemkit_test(graphicsray graphicsraytest)
