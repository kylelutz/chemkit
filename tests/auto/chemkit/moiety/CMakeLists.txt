qt4_wrap_cpp(MOC_SOURCES moietytest.h)
add_executable(moietytest moietytest.cpp ${MOC_SOURCES})
target_link_libraries(moietytest chemkit ${QT_LIBRARIES})
add_chemkit_test(chemkit.Moiety moietytest)
