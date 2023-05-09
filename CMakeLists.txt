cmake_minimum_required(VERSION 3.14)
project(hw3 LANGUAGES CXX)

find_package(UPCXX REQUIRED)

# Build the kmer_hash executable
add_executable(main_dist_pcsr main_dist_pcsr.cpp)
target_link_libraries(main_dist_pcsr PRIVATE UPCXX::upcxx)
# target_compile_definitions(main_dist_pcsr PRIVATE "KMER_LEN=19")

# Copy the job scripts
configure_file(job-perlmutter-starter job-perlmutter-starter COPYONLY)

