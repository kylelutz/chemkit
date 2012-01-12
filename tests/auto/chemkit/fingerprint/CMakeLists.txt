qt4_wrap_cpp(MOC_SOURCES fingerprinttest.h)
add_executable(fingerprinttest fingerprinttest.cpp ${MOC_SOURCES})
target_link_libraries(fingerprinttest chemkit ${QT_LIBRARIES})
add_chemkit_test(chemkit.Fingerprint fingerprinttest)
