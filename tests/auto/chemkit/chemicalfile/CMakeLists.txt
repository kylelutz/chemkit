qt4_wrap_cpp(MOC_SOURCES chemicalfiletest.h)
add_executable(chemicalfiletest chemicalfiletest.cpp ${MOC_SOURCES})
target_link_libraries(chemicalfiletest chemkit ${QT_LIBRARIES})
add_chemkit_test(chemicalfile chemicalfiletest)
