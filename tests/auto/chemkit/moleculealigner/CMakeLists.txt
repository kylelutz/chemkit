qt4_wrap_cpp(MOC_SOURCES moleculealignertest.h)
add_executable(moleculealignertest moleculealignertest.cpp ${MOC_SOURCES})
target_link_libraries(moleculealignertest chemkit chemkit-io ${QT_LIBRARIES})
add_chemkit_test(moleculealigner moleculealignertest)
