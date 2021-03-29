cmake_minimum_required(VERSION 2.6)
project(parametrization)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
set(EXECUTABLE_OUTPUT_PATH  ${CMAKE_CURRENT_BINARY_DIR}/bin)

### libIGL options: choose between header only and compiled static library
option(LIBIGL_USE_STATIC_LIBRARY     "Use libigl as static library" OFF)
# option(LIBIGL_WITH_OPENGL            "Use OpenGL"         ON)
# option(LIBIGL_WITH_OPENGL_GLFW       "Use GLFW"           ON)

find_package(LIBIGL REQUIRED)
find_package(VTK REQUIRED)


# Add your project files
include_directories("include/")
file(GLOB SRCFILES src/*.cpp)


add_executable(parametrization "main.cpp" ${SRCFILES})
target_link_libraries(parametrization                     
                    igl::core 
                    igl::opengl 
                    igl::opengl_glfw
                    ${VTK_LIBRARIES}
                    )
