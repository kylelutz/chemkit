qt4_wrap_cpp(MOC_SOURCES smilesinchitest.h)
add_executable(smilesinchitest smilesinchitest.cpp ${MOC_SOURCES})
target_link_libraries(smilesinchitest chemkit ${QT_LIBRARIES})
add_chemkit_test(conversion.SmilesInchi smilesinchitest)
