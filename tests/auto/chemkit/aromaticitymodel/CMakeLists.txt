qt4_wrap_cpp(MOC_SOURCES aromaticitymodeltest.h)
add_executable(aromaticitymodeltest aromaticitymodeltest.cpp mockaromaticitymodel.cpp ${MOC_SOURCES})
target_link_libraries(aromaticitymodeltest chemkit ${QT_LIBRARIES})
add_chemkit_test(chemkit.AromaticityModel aromaticitymodeltest)
