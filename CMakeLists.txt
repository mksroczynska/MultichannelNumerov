cmake_minimum_required(VERSION 3.5)

project(numerov)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -fopenmp -O3 ")

#set(SOURCE_FILES main.cpp Parameters.cpp Parameters.h ConstantGridSolver.cpp ConstantGridSolver.h NonconstantGridSolver.cpp NonconstantGridSolver.h Solver.cpp Solver.h)
set(SOURCE_FILES main.cpp Parameters.cpp Parameters.h ConstantGridSolver.cpp ConstantGridSolver.h NonconstantGridSolver.cpp NonconstantGridSolver.h)
add_executable(numerov ${SOURCE_FILES})
set(LIBS -larmadillo)
target_link_libraries(${PROJECT_NAME} ${LIBS} )

# Locate GTest
find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})




# Link run_tests with what we want to test and the GTest and pthread library
add_executable(run_tests Tests/tests.cpp Parameters.cpp Parameters.h ConstantGridSolver.cpp ConstantGridSolver.h Tests/ParametersTests.h Tests/ConstantGridSolverTests.h NonconstantGridSolverTests.h)
target_link_libraries(run_tests ${GTEST_LIBRARIES} ${LIBS} pthread)