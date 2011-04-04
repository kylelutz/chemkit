qt4_wrap_cpp(MOC_SOURCES ufftest.h)
add_executable(ufftest ufftest.cpp ${MOC_SOURCES})
target_link_libraries(ufftest chemkit ${QT_LIBRARIES})
add_chemkit_test(uff ufftest)
