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

# Headless, no SFML: runs without a display.
option(BUILD_TESTS "Build the host-only unit tests" ON)
if(BUILD_TESTS)
    add_executable(boundary_tests
            tests/test_boundary.cpp
            src/mesh.cpp
            src/solid.cpp)
    target_include_directories(boundary_tests PRIVATE
            include/
            /usr/include/suitesparse)
    target_link_libraries(boundary_tests PRIVATE
            Eigen3::Eigen
            OpenMP::OpenMP_CXX
            umfpack)
    target_compile_definitions(boundary_tests PRIVATE
            $<$<BOOL:${SHOULD_PROFILE}>:SHOULD_PROFILE>)
endif()

# Opt-in: needs python3 and resources/DXFs/cell2.dxf.
option(BUILD_DXF_DIAGNOSTICS "Build the headless real-DXF diagnostics" OFF)
if(BUILD_DXF_DIAGNOSTICS)
    add_executable(dxf_diagnostics
            tests/dxf_diagnostics.cpp
            src/mesh.cpp
            src/solid.cpp
            src/loader.cpp)
    target_include_directories(dxf_diagnostics PRIVATE
            include/
            /usr/include/suitesparse)
    target_compile_definitions(dxf_diagnostics PRIVATE
            SUBROUTINES_PATH="${CMAKE_SOURCE_DIR}/src/subroutines"
            DXFS_DIRECTORY_PATH="${CMAKE_SOURCE_DIR}/resources/DXFs")
    target_link_libraries(dxf_diagnostics PRIVATE
            Eigen3::Eigen
            OpenMP::OpenMP_CXX
            umfpack)
endif()
