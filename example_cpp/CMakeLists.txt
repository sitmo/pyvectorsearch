cmake_minimum_required(VERSION 3.2)
set (CMAKE_CXX_STANDARD 11)

# Set the project name
project(vectorsearch_tests)

# Create a sources variable with a link to all cpp files to compile
set(SOURCE_FILES main.cpp)

# Add an executable with the above sources
add_executable(vectorsearch_tests ${SOURCE_FILES})

include_directories(${PROJECT_SOURCE_DIR}/include)

# Set the directories that should be included in the build command for this target
# when running g++ these will be included as -I/directory/path/
target_include_directories(vectorsearch_tests PRIVATE ${PROJECT_SOURCE_DIR}/include)
install(TARGETS vectorsearch_tests DESTINATION bin)
