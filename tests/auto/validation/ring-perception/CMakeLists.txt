qt4_wrap_cpp(MOC_SOURCES ringperceptiontest.h)
add_executable(ringperceptiontest ringperceptiontest.cpp ${MOC_SOURCES})
target_link_libraries(ringperceptiontest chemkit ${QT_LIBRARIES})
add_chemkit_test(validation.RingPerception ringperceptiontest)
