find_package(Chemkit)
include_directories(${CHEMKIT_INCLUDE_DIRS})

if(NOT ${CHEMKIT_WITH_WIDGETS})
  return()
endif()

find_package(Qt4 4.6 COMPONENTS QtCore QtGui REQUIRED)
include(${QT_USE_FILE})

set(HEADERS
  atomcolormap.h
  periodictabledialog.h
  periodictablewidget.h
  widgets.h
)

set(SOURCES
  atomcolormap.cpp
  periodictabledialog.cpp
  periodictablewidget.cpp
  widgets.cpp
)

set(MOC_HEADERS
  periodictabledialog.h
  periodictablewidget.h
)

set(UI_FORMS
  periodictabledialog.ui
  periodictablewidget.ui
)

qt4_wrap_ui(UI_SOURCES ${UI_FORMS})
qt4_wrap_cpp(MOC_SOURCES ${MOC_HEADERS})

add_definitions(
  -DCHEMKIT_WIDGETS_LIBRARY
)

include_directories(${CMAKE_CURRENT_BINARY_DIR})

add_chemkit_library(chemkit-widgets ${SOURCES} ${UI_SOURCES} ${MOC_SOURCES})
target_link_libraries(chemkit-widgets ${CHEMKIT_LIBRARIES} ${QT_LIBRARIES})

# install header files
install(FILES ${HEADERS} DESTINATION include/chemkit/)
