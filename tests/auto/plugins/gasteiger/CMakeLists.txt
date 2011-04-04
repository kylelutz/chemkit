qt4_wrap_cpp(MOC_SOURCES gasteigertest.h)
add_executable(gasteigertest gasteigertest.cpp ${MOC_SOURCES})
target_link_libraries(gasteigertest chemkit ${QT_LIBRARIES})
add_chemkit_test(gasteiger gasteigertest)
