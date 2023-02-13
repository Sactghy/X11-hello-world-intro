cmake_minimum_required(VERSION 3.5)

project(Hello-c LANGUAGES C)

add_executable(Hello-c main.c)
target_link_libraries(Hello-c m)
target_link_libraries(Hello-c X11)
