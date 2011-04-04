qt4_wrap_cpp(MOC_SOURCES atommappingtest.h)
add_executable(atommappingtest atommappingtest.cpp ${MOC_SOURCES})
target_link_libraries(atommappingtest chemkit ${QT_LIBRARIES})
add_chemkit_test(atommapping atommappingtest)
