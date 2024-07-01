cmake_minimum_required(VERSION 3.22)
project(A_star_for_flipdistance)

set(CMAKE_CXX_STANDARD 17)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR})
find_package(GUROBI REQUIRED)

find_package(CGAL REQUIRED)
include(${CGAL_USE_FILE})


#set(CMAKE_CXX_FLAGS "-Wall -Wextra")
#set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

add_executable(A_star_for_flipdistance_experiments
        src/main.cpp
        src/Geometry.cpp
        src/ListHandler.cpp
        src/GeometricObjects.cpp
        src/TriangulationHandler.cpp
        src/StaticTriangulation.cpp
        src/AStarFlipDistance.cpp
        src/HungarianPrimDual.cpp
        src/HeuristicDistanceCalculater.cpp
        src/CombinatorialDelaunayTriangulation.cpp
        src/BiDirectionalBfs.cpp
        src/FastHankeHeuristic.cpp
        src/SimpleMaxHeap.cpp
        src/OuterBiConnectedComponentFinder.cpp
        src/Decomposer.cpp
        src/AStarDecompositionFlipDistance.cpp
        src/FastEppsteinHeuristic.cpp
)


#target_compile_definitions(A_star_for_flipdistance PUBLIC LARGEINTS=0)


target_include_directories(A_star_for_flipdistance_experiments
        PRIVATE
        ${PROJECT_SOURCE_DIR}/include
        ${GUROBI_INCLUDE_DIRS}
)

target_link_libraries(A_star_for_flipdistance_experiments
        optimized ${GUROBI_CXX_LIBRARY}
        debug ${GUROBI_CXX_DEBUG_LIBRARY}
        ${GUROBI_LIBRARY}
        CGAL::CGAL
)

add_executable(A_star_for_flipdistance_decomposition
        src/main.cpp
        src/Geometry.cpp
        src/ListHandler.cpp
        src/GeometricObjects.cpp
        src/TriangulationHandler.cpp
        src/StaticTriangulation.cpp
        src/AStarFlipDistance.cpp
        src/HungarianPrimDual.cpp
        src/HeuristicDistanceCalculater.cpp
        src/CombinatorialDelaunayTriangulation.cpp
        src/BiDirectionalBfs.cpp
        src/FastHankeHeuristic.cpp
        src/SimpleMaxHeap.cpp
        src/OuterBiConnectedComponentFinder.cpp
        src/Decomposer.cpp
        src/AStarDecompositionFlipDistance.cpp
        src/FastEppsteinHeuristic.cpp
)

target_compile_definitions(A_star_for_flipdistance_decomposition PUBLIC LARGEINTS=0)


target_include_directories(A_star_for_flipdistance_decomposition
        PRIVATE
        ${PROJECT_SOURCE_DIR}/include
        ${GUROBI_INCLUDE_DIRS}
)

target_link_libraries(A_star_for_flipdistance_decomposition
        optimized ${GUROBI_CXX_LIBRARY}
        debug ${GUROBI_CXX_DEBUG_LIBRARY}
        ${GUROBI_LIBRARY}
        CGAL::CGAL
)



add_executable(A_star_for_flipdistance
        src/main_with_path.cpp
        src/Geometry.cpp
        src/ListHandler.cpp
        src/GeometricObjects.cpp
        src/TriangulationHandler.cpp
        src/StaticTriangulation.cpp
        src/AStarFlipDistance.cpp
        src/HungarianPrimDual.cpp
        src/HeuristicDistanceCalculater.cpp
        src/CombinatorialDelaunayTriangulation.cpp
        src/BiDirectionalBfs.cpp
        src/FastHankeHeuristic.cpp
        src/SimpleMaxHeap.cpp
        src/OuterBiConnectedComponentFinder.cpp
        src/Decomposer.cpp
        src/AStarDecompositionFlipDistance.cpp
        src/FastEppsteinHeuristic.cpp
)


#target_compile_definitions(A_star_for_flipdistance_one_instance PUBLIC LARGEINTS=1)


target_include_directories(A_star_for_flipdistance
        PRIVATE
        ${PROJECT_SOURCE_DIR}/include
        ${GUROBI_INCLUDE_DIRS}
)

target_link_libraries(A_star_for_flipdistance
        optimized ${GUROBI_CXX_LIBRARY}
        debug ${GUROBI_CXX_DEBUG_LIBRARY}
        ${GUROBI_LIBRARY}
        CGAL::CGAL
)