cmake_minimum_required(VERSION 3.21)
project(snake_tracer CXX)
set(CMAKE_CXX_STANDARD 17)

find_package(OpenBLAS REQUIRED)

add_executable(matrix_power
    matrix_power/main.cpp
    matrix_power/matrix.cpp
)

target_include_directories(matrix_power PUBLIC 
    ${CMAKE_CURRENT_LIST_DIR}
)
target_link_libraries(matrix_power PUBLIC
    OpenBLAS::OpenBLAS
)