qt4_wrap_cpp(MOC_SOURCES hydrogenbondstest.h)
add_executable(hydrogenbondstest hydrogenbondstest.cpp ${MOC_SOURCES})
target_link_libraries(hydrogenbondstest chemkit ${QT_LIBRARIES})
add_chemkit_test(plugins.HydrogenBonds hydrogenbondstest)
