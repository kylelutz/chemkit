qt4_wrap_cpp(MOC_SOURCES graphicsscenetest.h)
add_executable(graphicsscenetest graphicsscenetest.cpp ${MOC_SOURCES})
target_link_libraries(graphicsscenetest chemkit chemkit-graphics ${QT_LIBRARIES})
add_chemkit_test(graphicsscene graphicsscenetest)
