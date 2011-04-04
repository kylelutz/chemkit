qt4_wrap_cpp(MOC_SOURCES alphashapetest.h)
add_executable(alphashapetest alphashapetest.cpp ${MOC_SOURCES})
target_link_libraries(alphashapetest chemkit ${QT_LIBRARIES})
add_chemkit_test(alphashape alphashapetest)
