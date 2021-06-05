cmake_minimum_required(VERSION 3.18cc)
project(parametrization)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
set(EXECUTABLE_OUTPUT_PATH  ${CMAKE_CURRENT_BINARY_DIR}/bin)

### libIGL options: choose between header only and compiled static library
# option(LIBIGL_USE_STATIC_LIBRARY     "Use libigl as static library" ON)


find_package(LIBIGL REQUIRED)
find_package(VTK REQUIRED)


# Add your project files
include_directories("include/")
include_directories("libigl/include")
include_directories("libigl/external/eigen")

file(GLOB SRCFILES src/*.cpp)


add_executable(parametrization "main.cpp" ${SRCFILES})
target_link_libraries(parametrization                     
                    # igl::core 
                    ${VTK_LIBRARIES}
                    )
