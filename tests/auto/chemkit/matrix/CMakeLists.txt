qt4_wrap_cpp(MOC_SOURCES matrixtest.h)
add_executable(matrixtest matrixtest.cpp ${MOC_SOURCES})
target_link_libraries(matrixtest chemkit ${QT_LIBRARIES})
add_chemkit_test(matrix matrixtest)
