cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
#set(CMAKE_CXX_FLAGS "-Wall -pedantic -Wextra")

project(kinectv2_decoders)
find_package(HDF5 REQUIRED COMPONENTS C CXX HL)

include_directories(include ${KINECTV2_DECODERS_INCLUDE_DIRS} ${HDF5_INCLUDE_DIRS})

set(kinectv2_decoders_src run_decoders.cpp
                        src/libfreenect2_processor.cpp
                        src/kde_processor.cpp
                        src/device_parameters.cpp
                        src/read_file_handler.cpp
                        src/setup_parser.cpp)

add_executable(kinectv2_decoders ${kinectv2_decoders_src})
target_link_libraries(kinectv2_decoders ${HDF5_LIBRARIES} ${HDF5_HL_LIBRARIES})

