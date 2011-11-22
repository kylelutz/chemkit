qt4_wrap_cpp(MOC_SOURCES massdescriptorstest.h)
add_executable(massdescriptorstest massdescriptorstest.cpp ${MOC_SOURCES})
target_link_libraries(massdescriptorstest chemkit ${QT_LIBRARIES})
add_chemkit_test(plugins.MassDescriptors massdescriptorstest)
