cmake_minimum_required(VERSION 3.16)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
project(runThis)

set(CMAKE_CXX_STANDARD 17)

file(GLOB_RECURSE SOURCES "src/*.cpp")

add_executable(${PROJECT_NAME} ${SOURCES}
        include/constants.hpp
        include/mesh.hpp
        include/renderer.hpp
        include/solid.hpp
        include/loader.hpp
        include/utils.hpp
        src/loader.cpp
)

# find_package(OpenGL REQUIRED)
# find_package(glfw3 3.3 REQUIRED)
target_compile_definitions(${PROJECT_NAME} PRIVATE
        PARTICLE_TEXTURE="${CMAKE_SOURCE_DIR}/resources/textures/circle.png"

)
target_include_directories(${PROJECT_NAME} PRIVATE
        include/
)

# find_package(glm REQUIRED)
# find_package(assimp REQUIRED)

# target_link_libraries(${PROJECT_NAME} PRIVATE
#     OpenGL::GL
#     glfw
#     glm::glm
#     assimp::assimp
# )