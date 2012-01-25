qt4_wrap_cpp(MOC_SOURCES atomcolormaptest.h)
add_executable(atomcolormaptest atomcolormaptest.cpp ${MOC_SOURCES})
target_link_libraries(atomcolormaptest chemkit chemkit-gui ${QT_LIBRARIES})
add_chemkit_test(gui.AtomColorMap atomcolormaptest)
