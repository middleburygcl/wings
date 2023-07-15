# add all the source files
set(SOURCES
  field.cpp
  io.cpp
  log.cpp
  mesh.cpp
  shader.cpp
  util.cpp)

# create an object library for all the sources
add_library(wings_util_src OBJECT ${SOURCES})
add_definitions(-DWINGS_SOURCE_DIR="${PROJECT_SOURCE_DIR}")
target_compile_definitions(wings_util_src PUBLIC -DWINGS_SOURCE_DIR="${PROJECT_SOURCE_DIR}")

# create both shared and static libraries with same name
add_library(wings_util_shared SHARED $<TARGET_OBJECTS:wings_util_src>)
set_target_properties(wings_util_shared PROPERTIES OUTPUT_NAME wings_util CLEAN_DIRECT_OUTPUT 1)
target_link_libraries(wings_util_shared dl)
if (UNIX AND NOT APPLE)
  target_link_libraries(wings_util_shared rt dl)
endif()
add_library(wings_util_static STATIC $<TARGET_OBJECTS:wings_util_src>)
set_target_properties(wings_util_static PROPERTIES OUTPUT_NAME wings_util CLEAN_DIRECT_OUTPUT 1)

target_link_libraries(wings_util_shared ${WINGS_EXTERNAL_LIBRARIES})
target_link_libraries(wings_util_static ${WINGS_EXTERNAL_LIBRARIES})


