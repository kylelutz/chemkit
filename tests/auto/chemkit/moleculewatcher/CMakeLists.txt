qt4_wrap_cpp(MOC_SOURCES moleculewatchertest.h)
add_executable(moleculewatchertest moleculewatchertest.cpp ${MOC_SOURCES})
target_link_libraries(moleculewatchertest chemkit ${QT_LIBRARIES})
add_chemkit_test(chemkit.MoleculeWatcher moleculewatchertest)
