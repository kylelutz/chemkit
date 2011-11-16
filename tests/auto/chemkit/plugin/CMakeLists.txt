qt4_wrap_cpp(MOC_SOURCES plugintest.h)
add_executable(plugintest plugintest.cpp mockplugin.cpp ${MOC_SOURCES})
target_link_libraries(plugintest chemkit ${QT_LIBRARIES})
add_chemkit_test(chemkit.Plugin plugintest)
