qt4_wrap_cpp(MOC_SOURCES inchitest.h)
add_executable(inchitest inchitest.cpp ${MOC_SOURCES})
target_link_libraries(inchitest chemkit chemkit-io ${QT_LIBRARIES})
add_chemkit_test(inchi inchitest)
