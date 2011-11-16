qt4_wrap_cpp(MOC_SOURCES bondpredictortest.h)
add_executable(bondpredictortest bondpredictortest.cpp ${MOC_SOURCES})
target_link_libraries(bondpredictortest chemkit ${QT_LIBRARIES})
add_chemkit_test(chemkit.BondPredictor bondpredictortest)
