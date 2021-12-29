cmake_minimum_required(VERSION 3.18)
project(ca LANGUAGES CXX C)

set(CMAKE_CXX_STANDARD 20)
if (MSVC)
    add_compile_options(/W4 /WX)
else()
    add_compile_options(-Wall -Wextra -pedantic -Werror)
endif()

add_subdirectory(src)
add_subdirectory(tests)
