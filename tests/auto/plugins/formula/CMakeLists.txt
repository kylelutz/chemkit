qt4_wrap_cpp(MOC_SOURCES formulatest.h)
add_executable(formulatest formulatest.cpp ${MOC_SOURCES})
target_link_libraries(formulatest chemkit ${QT_LIBRARIES})
add_chemkit_test(formula formulatest)
