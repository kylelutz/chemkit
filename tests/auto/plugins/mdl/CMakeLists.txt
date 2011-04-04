qt4_wrap_cpp(MOC_SOURCES mdltest.h)
add_executable(mdltest mdltest.cpp ${MOC_SOURCES})
target_link_libraries(mdltest chemkit ${QT_LIBRARIES})
add_chemkit_test(mdl mdltest)
