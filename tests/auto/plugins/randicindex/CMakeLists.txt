qt4_wrap_cpp(MOC_SOURCES randicindextest.h)
add_executable(randicindextest randicindextest.cpp ${MOC_SOURCES})
target_link_libraries(randicindextest chemkit ${QT_LIBRARIES})
add_chemkit_test(plugins.RandicIndex randicindextest)
