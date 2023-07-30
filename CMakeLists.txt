cmake_minimum_required(VERSION 3.14)
project(hw3 LANGUAGES CXX)

find_package(UPCXX REQUIRED)

add_executable(main_dist_pcsr src/main/main_dist_pcsr.cpp)
target_link_libraries(main_dist_pcsr PRIVATE UPCXX::upcxx)
# target_compile_definitions(main_dist_pcsr PRIVATE "KMER_LEN=19")

