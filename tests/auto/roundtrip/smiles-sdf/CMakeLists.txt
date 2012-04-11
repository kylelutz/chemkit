qt4_wrap_cpp(MOC_SOURCES smilessdftest.h)
add_executable(smilessdftest smilessdftest.cpp ${MOC_SOURCES})
target_link_libraries(smilessdftest chemkit chemkit-io ${QT_LIBRARIES})
add_chemkit_test(roundtrip.SmilesSdf smilessdftest)
