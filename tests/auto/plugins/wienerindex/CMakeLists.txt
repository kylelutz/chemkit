qt4_wrap_cpp(MOC_SOURCES wienerindextest.h)
add_executable(wienerindextest wienerindextest.cpp ${MOC_SOURCES})
target_link_libraries(wienerindextest chemkit ${QT_LIBRARIES})
add_chemkit_test(plugins.WienerIndex wienerindextest)
