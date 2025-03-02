cmake_minimum_required(VERSION 3.27)
project(flip)
set(CMAKE_CXX_STANDARD 23)
set(SOURCE_DIR src)
set(VENDOR_DIR vendor)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(EMSCRIPTEN)
  add_compile_definitions(PLATFORM_WEB)
  set(BUILD_PLATFORM "PLATFORM=Web")
else()
  add_compile_definitions(PLATFORM_DESKTOP)
  set(BUILD_PLATFORM "PLATFORM=Desktop")
endif()

add_subdirectory("vendor/raylib")

file(GLOB_RECURSE SOURCE_FILES
  "${SOURCE_DIR}/**.cpp"
  "${SOURCE_DIR}/**.h"
)

set(INCLUDE_DIRECTORIES
  "${SOURCE_DIR}"
  "${VENDOR_DIR}/raylib/src"
)

add_executable(flip "${SOURCE_FILES}")
target_include_directories(flip PUBLIC ${INCLUDE_DIRECTORIES})
target_link_libraries(flip PUBLIC raylib)

if(EMSCRIPTEN)
  set_target_properties(flip PROPERTIES
    SUFFIX ".js"
    LINK_FLAGS "-s USE_GLFW=3 -s INITIAL_MEMORY=16777216 -s ASYNCIFY=1 -s ASYNCIFY_IMPORTS=['emscripten_sleep'] -s EXPORT_NAME='flip' -s EXPORTED_RUNTIME_METHODS=['ccall','cwrap'] --shell-file ${CMAKE_SOURCE_DIR}/index.html --preload-file ${CMAKE_SOURCE_DIR}/assets@assets --preload-file ${CMAKE_SOURCE_DIR}/shaders@shaders"
  )
  
  # Make sure these files get copied to the build directory
  add_custom_command(TARGET flip POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy 
      ${CMAKE_SOURCE_DIR}/index.html
      ${CMAKE_BINARY_DIR}/index.html
  )
  
  add_custom_command(TARGET flip POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
      ${CMAKE_SOURCE_DIR}/assets
      ${CMAKE_BINARY_DIR}/assets
    COMMAND ${CMAKE_COMMAND} -E copy_directory
      ${CMAKE_SOURCE_DIR}/shaders
      ${CMAKE_BINARY_DIR}/shaders
  )
endif()