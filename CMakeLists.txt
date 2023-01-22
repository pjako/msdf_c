cmake_minimum_required(VERSION 3.19)
project(engine)

add_definitions(-DSAMPLE_ROOT="${PROJECT_SOURCE_DIR}/sample")

set(CMAKE_C_STANDARD 11)

add_executable(msdf_sample sample/sample.c)