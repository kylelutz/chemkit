qt4_wrap_cpp(MOC_SOURCES substructurequerytest.h)
add_executable(substructurequerytest substructurequerytest.cpp ${MOC_SOURCES})
target_link_libraries(substructurequerytest chemkit ${QT_LIBRARIES})
add_chemkit_test(chemkit.SubstructureQuery substructurequerytest)
