cmake_minimum_required(VERSION 3.16)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
project(runThis)

set(CMAKE_CXX_STANDARD 17)

file(GLOB_RECURSE SOURCES "src/*.cpp")
set(IMGUI_DIR ${CMAKE_SOURCE_DIR}/external/imgui)
set(IMGUI_SFML_DIR ${CMAKE_SOURCE_DIR}/external/imgui-sfml)

file(GLOB IMGUI_SOURCES ${IMGUI_DIR}/*.cpp)

add_executable(${PROJECT_NAME} ${SOURCES}
        include/constants.hpp
        include/mesh.hpp
        include/renderer.hpp
        include/solid.hpp
        include/loader.hpp
        include/utils.hpp
        src/loader.cpp
        ${IMGUI_SOURCES}
        ${IMGUI_SFML_DIR}/imgui-SFML.cpp
        include/menu.hpp
        src/menu.cpp
        include/flag.hpp
)

# find_package(OpenGL REQUIRED)
# find_package(glfw3 3.3 REQUIRED)
target_compile_definitions(${PROJECT_NAME} PRIVATE
        PARTICLE_TEXTURE_PATH="${CMAKE_SOURCE_DIR}/resources/textures/circle.png"
        DXFS_DIRECTORY_PATH="${CMAKE_SOURCE_DIR}/resources/DXFs"
        SUBROUTINES_PATH="${CMAKE_SOURCE_DIR}/src/subroutines"

)
target_include_directories(${PROJECT_NAME} PRIVATE
        include/
        ${IMGUI_DIR}
        ${IMGUI_SFML_DIR}
)

# find_package(glm REQUIRED)
# find_package(assimp REQUIRED)
find_package(SFML 3 COMPONENTS Graphics REQUIRED)

# target_link_libraries(${PROJECT_NAME} PRIVATE
#     OpenGL::GL
#     glfw
#     glm::glm
#     assimp::assimp
# )
target_link_libraries(${PROJECT_NAME} PRIVATE
        SFML::Graphics
        SFML::Window
        SFML::System
)

find_package(OpenGL REQUIRED)
target_link_libraries(${PROJECT_NAME} PRIVATE OpenGL::GL)

