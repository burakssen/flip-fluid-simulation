cmake_minimum_required(VERSION 3.20)  
project(flip)


set(CMAKE_CXX_STANDARD 23)

set(OMP_DIR "/opt/homebrew/Cellar/libomp/19.1.7")

# Add Raylib
add_subdirectory(vendor/raylib)

# Add source files
file(GLOB_RECURSE SOURCES "src/*.cpp")
add_executable(${PROJECT_NAME} ${SOURCES})

# Configure include directories
target_include_directories(${PROJECT_NAME} 
    PRIVATE 
        vendor/raylib/src
        ${OMP_DIR}/include
)

# Ensure correct linking with ANGLE
target_link_libraries(${PROJECT_NAME} 
    PRIVATE 
        raylib
        ${OMP_DIR}/lib/libomp.a
)
