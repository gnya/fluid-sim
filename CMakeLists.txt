cmake_minimum_required(VERSION 3.10)

project(fluid_sim)

enable_language(CXX)
set(CMAKE_CXX_STANDARD 14)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -mavx -lpng")

set(CMAKE_THREAD_PREFER_PTHREAD TRUE)
set(THREADS_PREFER_PTHREAD_FLAG TRUE)
set(OpenGL_GL_PREFERENCE GLVND)

find_package(Threads REQUIRED)
find_package(GLUT REQUIRED)
find_package(OpenGL REQUIRED)
find_package(PNG REQUIRED)
find_package(Boost COMPONENTS program_options REQUIRED)

include_directories(${OPENGL_INCLUDE_DIRS} ${GLUT_INCLUDE_DIRS} ${PNG_INCLUDE_DIR}  ${BOOST_INCLUDE_DIR})

add_executable(fluid_sim main.cpp lib/fluid.h lib/map_fluid.h lib/noise.h lib/ink_fluid.h)
target_link_libraries(fluid_sim Threads::Threads GLUT::GLUT OpenGL::GL PNG::PNG ${Boost_LIBRARIES})