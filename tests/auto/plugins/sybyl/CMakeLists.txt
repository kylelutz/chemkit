qt4_wrap_cpp(MOC_SOURCES sybyltest.h)
add_executable(sybyltest sybyltest.cpp ${MOC_SOURCES})
target_link_libraries(sybyltest chemkit ${QT_LIBRARIES})
add_chemkit_test(sybyl sybyltest)
