set(QT_USE_QTXML TRUE)
include(${QT_USE_FILE})

qt4_wrap_cpp(MOC_SOURCES mmfftest.h)
add_executable(mmfftest mmfftest.cpp ${MOC_SOURCES})
target_link_libraries(mmfftest chemkit chemkit-io ${QT_LIBRARIES})
add_chemkit_test(mmff mmfftest)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/mmff94.expected DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
