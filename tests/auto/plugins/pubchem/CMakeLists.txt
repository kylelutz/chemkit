qt4_wrap_cpp(MOC_SOURCES pubchemtest.h)
add_executable(pubchemtest pubchemtest.cpp ${MOC_SOURCES})
target_link_libraries(pubchemtest chemkit ${QT_LIBRARIES})
add_chemkit_test(plugins.PubChem pubchemtest)
