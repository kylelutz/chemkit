qt4_wrap_cpp(MOC_SOURCES txyztest.h)
add_executable(txyztest txyztest.cpp ${MOC_SOURCES})
target_link_libraries(txyztest chemkit ${QT_LIBRARIES})
add_chemkit_test(txyz txyztest)
