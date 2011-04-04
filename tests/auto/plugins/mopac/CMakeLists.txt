qt4_wrap_cpp(MOC_SOURCES mopactest.h)
add_executable(mopactest mopactest.cpp ${MOC_SOURCES})
target_link_libraries(mopactest chemkit ${QT_LIBRARIES})
add_chemkit_test(mopac mopactest)
