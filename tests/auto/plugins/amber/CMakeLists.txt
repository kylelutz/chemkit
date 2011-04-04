qt4_wrap_cpp(MOC_SOURCES ambertest.h)
add_executable(ambertest ambertest.cpp ${MOC_SOURCES})
target_link_libraries(ambertest chemkit ${QT_LIBRARIES})
add_chemkit_test(amber ambertest)
