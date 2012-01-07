qt4_wrap_cpp(MOC_SOURCES sdfmol2test.h)
add_executable(sdfmol2test sdfmol2test.cpp ${MOC_SOURCES})
target_link_libraries(sdfmol2test chemkit chemkit-io ${QT_LIBRARIES})
add_chemkit_test(conversion.SdfMol2 sdfmol2test)
