qt4_wrap_cpp(MOC_SOURCES apoltest.h)
add_executable(apoltest apoltest.cpp ${MOC_SOURCES})
target_link_libraries(apoltest chemkit ${QT_LIBRARIES})
add_chemkit_test(plugins.Apol apoltest)
