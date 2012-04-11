find_package(Chemkit REQUIRED)
include_directories(${CHEMKIT_INCLUDE_DIRS})

set(SOURCES
  formulalineformat.cpp
  formulaplugin.cpp
  htmlformulalineformat.cpp
  spacedformulalineformat.cpp
)

add_chemkit_plugin(formula ${SOURCES})
target_link_libraries(formula ${CHEMKIT_LIBRARIES})
