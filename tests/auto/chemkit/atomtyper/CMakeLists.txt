qt4_wrap_cpp(MOC_SOURCES atomtypertest.h)
add_executable(atomtypertest atomtypertest.cpp mockatomtyper.cpp ${MOC_SOURCES})
target_link_libraries(atomtypertest chemkit ${QT_LIBRARIES})
add_chemkit_test(chemkit.AtomTyper atomtypertest)
