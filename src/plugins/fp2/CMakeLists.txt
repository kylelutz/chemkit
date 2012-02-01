find_package(Chemkit REQUIRED)
include_directories(${CHEMKIT_INCLUDE_DIRS})

set(SOURCES
  fp2fingerprint.cpp
  fp2plugin.cpp
)

add_chemkit_plugin(fp2 ${SOURCES})
target_link_libraries(fp2 ${CHEMKIT_LIBRARIES})
