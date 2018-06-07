cmake_minimum_required(VERSION 3.5)

# Project name shown when building
set(PROJECT_NAME Assignment_5_Project)

# Directories specified here, no hardcoding
set(EXTERN_DIR ${CMAKE_CURRENT_SOURCE_DIR}/extern)
set(DATA_DIR ${CMAKE_CURRENT_SOURCE_DIR}/data)
set(OUTPUT_DIR ${CMAKE_CURRENT_SOURCE_DIR})
set(INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include)
set(SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)

# Create project
project(${PROJECT_NAME})

# Compiler and IDE configuration
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set_property(GLOBAL PROPERTY USE_FOLDERS ON)
add_definitions(-DDATA_DIR="${DATA_DIR}")
add_definitions(-DOUTPUT_DIR="${OUTPUT_DIR}")

# Required packages
find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)

# Glob all available files in include and src dirs. Kind of a bad practice...
set(INCLUDE_T ${INCLUDE_DIR}/*.hpp ${INCLUDE_DIR}/*.h)
set(SRC_T ${SRC_DIR}/*.cpp ${SRC_DIR}/*.c)
file(GLOB INCLUDES ${INCLUDE_T})
file(GLOB SRCS ${SRC_T})

# Create executable and link everything to it
add_executable(${PROJECT_NAME} ${SRCS} ${INCLUDES})
target_link_libraries(${PROJECT_NAME} ${OPENGL_LIBRARIES} ${GLUT_LIBRARIES})
target_include_directories(${PROJECT_NAME} PUBLIC 
    ${EXTERN_DIR}
    ${OPENGL_INCLUDE_DIR} 
    ${GLUT_INCLUDE_DIR} 
    ${INCLUDE_DIR})