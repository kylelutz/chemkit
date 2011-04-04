qt4_wrap_cpp(MOC_SOURCES scalarfieldtest.h)
add_executable(scalarfieldtest scalarfieldtest.cpp ${MOC_SOURCES})
target_link_libraries(scalarfieldtest chemkit ${QT_LIBRARIES})
add_chemkit_test(scalarfield scalarfieldtest)
