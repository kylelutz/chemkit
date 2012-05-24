qt4_wrap_cpp(MOC_SOURCES mannholdtest.h)
add_executable(mannholdtest mannholdtest.cpp ${MOC_SOURCES})
target_link_libraries(mannholdtest chemkit ${QT_LIBRARIES})
add_chemkit_test(plugins.Mannhold mannholdtest)
