qt4_wrap_cpp(MOC_SOURCES fhztest.h)
add_executable(fhztest fhztest.cpp ${MOC_SOURCES})
target_link_libraries(fhztest chemkit ${QT_LIBRARIES})
add_chemkit_test(fhz fhztest)
