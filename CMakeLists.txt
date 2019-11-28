cmake_minimum_required(VERSION 3.10)

project(fluid_sim)

enable_language(CXX)
set(CMAKE_CXX_STANDARD 17)

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=native")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -march=native")

set(CMAKE_THREAD_PREFER_PTHREAD TRUE)
set(THREADS_PREFER_PTHREAD_FLAG TRUE)
set(OpenGL_GL_PREFERENCE GLVND)

find_package(Threads REQUIRED)
find_package(PNG REQUIRED)
find_package(GLUT REQUIRED)
find_package(OpenGL REQUIRED)
find_package(Boost COMPONENTS program_options REQUIRED)

include_directories(${OPENGL_INCLUDE_DIRS})
include_directories(${GLUT_INCLUDE_DIRS})
include_directories(${PNG_INCLUDE_DIR})
include_directories(${BOOST_INCLUDE_DIR})
include_directories("${PROJECT_SOURCE_DIR}/lib/")

add_executable(fluid_sim main.cpp)
target_link_libraries(fluid_sim Threads::Threads PNG::PNG GLUT::GLUT OpenGL::GL ${Boost_LIBRARIES})