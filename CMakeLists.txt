cmake_minimum_required(VERSION 3.5)

project(EqnConc LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

include_directories(${PROJECT_SOURCE_DIR} ${PROJECT_BINARY_DIR} )


if(NOT EXISTS ${PROJECT_BINARY_DIR}/json.hpp)
  file(DOWNLOAD https://github.com/nlohmann/json/releases/download/v3.9.1/json.hpp ${PROJECT_BINARY_DIR}/json.hpp
       TIMEOUT 60  # seconds
       EXPECTED_HASH SHA256=7804b38146921d03374549c9e2a5e3acda097814c43caf2b96a0278e58df26e0
       TLS_VERIFY ON)
endif()

add_executable(EqnConc src/main.cpp)
target_link_libraries(EqnConc )
