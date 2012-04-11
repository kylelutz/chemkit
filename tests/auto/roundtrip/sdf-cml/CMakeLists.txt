qt4_wrap_cpp(MOC_SOURCES sdfcmltest.h)
add_executable(sdfcmltest sdfcmltest.cpp ${MOC_SOURCES})
target_link_libraries(sdfcmltest chemkit chemkit-io ${QT_LIBRARIES})
add_chemkit_test(roundtrip.SdfCml sdfcmltest)
