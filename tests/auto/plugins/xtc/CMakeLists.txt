qt4_wrap_cpp(MOC_SOURCES xtctest.h)
add_executable(xtctest xtctest.cpp ${MOC_SOURCES})
target_link_libraries(xtctest chemkit ${QT_LIBRARIES})
add_chemkit_test(xtc xtctest)
