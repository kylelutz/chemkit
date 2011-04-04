qt4_wrap_cpp(MOC_SOURCES rotatablebondstest.h)
add_executable(rotatablebondstest rotatablebondstest.cpp ${MOC_SOURCES})
target_link_libraries(rotatablebondstest chemkit ${QT_LIBRARIES})
add_chemkit_test(rotatablebonds rotatablebondstest)
