cmake_minimum_required(VERSION 3.10)
project(cagenet_preprocess)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Finde systemweit installiertes Eigen
find_package(Eigen3 REQUIRED)

# Füge deine Quellcodes hinzu
add_executable(mvc mvc.cpp)

# Linke Eigen (header-only)
target_link_libraries(mvc Eigen3::Eigen)

# Include-Pfade automatisch setzen
target_include_directories(mvc PRIVATE ${EIGEN3_INCLUDE_DIR})
