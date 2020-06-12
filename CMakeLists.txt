cmake_minimum_required(VERSION 3.15)
project(two_dimensional_heat_equation)

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_BUILD_TYPE "Debug")

if(CMAKE_BUILD_TYPE MATCHES Debug)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g3 -Wall -Wno-long-long -pedantic")
else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g0 -march=native -Os")
endif()

find_package(MPI REQUIRED)

include_directories(${MPI_INCLUDE_PATH})

aux_source_directory(src SRC)

add_executable(${PROJECT_NAME} ${SRC})

if(MPI_COMPILE_FLAGS)
    set_target_properties(${PROJECT_NAME} PROPERTIES
            COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
    set_target_properties(${PROJECT_NAME} PROPERTIES
            LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()

target_link_libraries(${PROJECT_NAME} ${MPI_C_LIBRARIES} ${MPI_CXX_LIBRARIES})

add_subdirectory(test)