cmake_minimum_required(VERSION 3.14)
project(hw3 LANGUAGES CXX)

find_package(UPCXX REQUIRED)

# Group number
set(GROUP_NAME "None" CACHE STRING "Your group name as it appears on bCourses (no spaces)")

# Packaging system
if (NOT ${GROUP_NAME} STREQUAL None)
    set(CPACK_GENERATOR TGZ)
    set(CPACK_PACKAGE_FILE_NAME "cs267${GROUP_NAME}_hw3")
    install(FILES kmer_hash.cpp hash_map.hpp DESTINATION .)
    install(FILES ${CPACK_PACKAGE_FILE_NAME}.pdf DESTINATION .)
    include(CPack)
endif ()

# Build the kmer_hash executable
add_executable(main_dist_pcsr main_dist_pcsr.cpp)
target_link_libraries(main_dist_pcsr PRIVATE UPCXX::upcxx)
target_compile_definitions(main_dist_pcsr PRIVATE "KMER_LEN=19")

# Copy the job scripts
configure_file(job-perlmutter-starter job-perlmutter-starter COPYONLY)

