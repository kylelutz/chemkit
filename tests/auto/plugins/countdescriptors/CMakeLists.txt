qt4_wrap_cpp(MOC_SOURCES countdescriptorstest.h)
add_executable(countdescriptorstest countdescriptorstest.cpp ${MOC_SOURCES})
target_link_libraries(countdescriptorstest chemkit ${QT_LIBRARIES})
add_chemkit_test(plugins.CountDescriptors countdescriptorstest)
