cmake_minimum_required(VERSION 3.16)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
project(runThis)

set(CMAKE_CXX_STANDARD 20)
option(SHOULD_PROFILE "Enable runtime profiling output" OFF)

file(GLOB_RECURSE SOURCES "src/*.cpp")
set(IMGUI_DIR ${CMAKE_SOURCE_DIR}/external/imgui)
set(IMGUI_SFML_DIR ${CMAKE_SOURCE_DIR}/external/imgui-sfml)
set(IMPLOT_DIR ${CMAKE_SOURCE_DIR}/external/implot)

file(GLOB IMGUI_SOURCES ${IMGUI_DIR}/*.cpp)
file(GLOB IMPLOT_SOURCES ${IMPLOT_DIR}/*.cpp)

add_executable(${PROJECT_NAME} ${SOURCES}
        include/constants.hpp
        include/mesh.hpp
        include/renderer.hpp
        include/solid.hpp
        include/loader.hpp
        include/utils.hpp
        src/loader.cpp
        ${IMGUI_SOURCES}
        ${IMPLOT_SOURCES}
        ${IMGUI_SFML_DIR}/imgui-SFML.cpp
        include/menu.hpp
        src/menu.cpp
        include/flag.hpp
        include/profiling.hpp
        include/vec.hpp
        include/domains.hpp
)
# find_package(OpenGL REQUIRED)
# find_package(glfw3 3.3 REQUIRED)
target_compile_definitions(${PROJECT_NAME} PRIVATE
        PARTICLE_TEXTURE_PATH="${CMAKE_SOURCE_DIR}/resources/textures/circle.png"
        DXFS_DIRECTORY_PATH="${CMAKE_SOURCE_DIR}/resources/DXFs"
        SUBROUTINES_PATH="${CMAKE_SOURCE_DIR}/src/subroutines"
        OUTPUT_VIDEOS_DIR="${CMAKE_SOURCE_DIR}/output/videos"
        OUTPUT_IMAGES_DIR="${CMAKE_SOURCE_DIR}/output/images"
        $<$<BOOL:${SHOULD_PROFILE}>:SHOULD_PROFILE>

)



find_package(SFML 3 COMPONENTS Graphics REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(OpenMP REQUIRED)

# find_package(PkgConfig REQUIRED)
# pkg_check_modules(UMFPACK REQUIRED umfpack)


target_include_directories(${PROJECT_NAME} PRIVATE
        include/
        ${IMGUI_DIR}
        ${IMPLOT_DIR}
        ${IMGUI_SFML_DIR}
        /usr/include/suitesparse
)

target_link_libraries(${PROJECT_NAME} PRIVATE
    SFML::Graphics
    SFML::Window
    SFML::System
    Eigen3::Eigen
    OpenMP::OpenMP_CXX
    OpenGL::GL

    umfpack
)

find_package(OpenGL REQUIRED)
target_link_libraries(${PROJECT_NAME} PRIVATE OpenGL::GL)

