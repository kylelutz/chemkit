find_package(Chemkit COMPONENTS io REQUIRED)
include_directories(${CHEMKIT_INCLUDE_DIRS})

set(SOURCES
  inchikeylineformat.cpp
  inchilineformat.cpp
  inchiplugin.cpp
)

aux_source_directory(../../3rdparty/inchi/ IUPAC_INCHI_SOURCES)

add_chemkit_plugin(inchi ${SOURCES} ${IUPAC_INCHI_SOURCES})
target_link_libraries(inchi ${CHEMKIT_LIBRARIES})
