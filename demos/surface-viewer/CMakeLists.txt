if(NOT ${CHEMKIT_WITH_IO} OR NOT ${CHEMKIT_WITH_GRAPHICS})
  return()
endif()

find_package(Chemkit COMPONENTS io graphics)
include_directories(${CHEMKIT_INCLUDE_DIRS})

find_package(Qt4 4.6 COMPONENTS QtCore QtGui QtOpenGL REQUIRED)
set(QT_USE_QTOPENGL TRUE)
include(${QT_USE_FILE})

include_directories(${CMAKE_CURRENT_BINARY_DIR})

set(SOURCES surfaceviewerdemo.cpp)
set(MOC_HEADERS surfaceviewerdemo.h)
set(UI_FORMS surfaceviewerdemo.ui)

qt4_wrap_ui(UI_SOURCES ${UI_FORMS})
qt4_wrap_cpp(MOC_SOURCES ${MOC_HEADERS})

add_executable(surface-viewer ${SOURCES} ${MOC_SOURCES} ${UI_SOURCES})
target_link_libraries(surface-viewer ${CHEMKIT_LIBRARIES} ${QT_LIBRARIES})
