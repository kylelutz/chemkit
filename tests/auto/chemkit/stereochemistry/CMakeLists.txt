qt4_wrap_cpp(MOC_SOURCES stereochemistrytest.h)
add_executable(stereochemistrytest stereochemistrytest.cpp ${MOC_SOURCES})
target_link_libraries(stereochemistrytest chemkit ${QT_LIBRARIES})
add_chemkit_test(stereochemistry stereochemistrytest)
