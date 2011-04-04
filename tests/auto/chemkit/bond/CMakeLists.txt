qt4_wrap_cpp(MOC_SOURCES bondtest.h)
add_executable(bondtest bondtest.cpp ${MOC_SOURCES})
target_link_libraries(bondtest chemkit ${QT_LIBRARIES})
add_chemkit_test(bond bondtest)
