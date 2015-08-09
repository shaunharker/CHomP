cmake_minimum_required (VERSION 2.6) 
project (CHomP) 

#########
# rpath #
#########

set(CMAKE_MACOSX_RPATH 1)
set(CMAKE_SKIP_BUILD_RPATH  FALSE)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE) 
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
if("${isSystemDir}" STREQUAL "-1")
   set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
endif("${isSystemDir}" STREQUAL "-1")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake-modules)

#########
# Flags #
#########

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -ftemplate-depth=2048")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 -fprofile-arcs -ftest-coverage")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -DNDEBUG" )

############
# Packages #
############

find_package(Boost 1.58.0 REQUIRED COMPONENTS chrono system thread serialization)
find_package(Threads)
find_package(X11 REQUIRED)
find_package(GMP REQUIRED)

#######################
# Include Directories #
#######################

include_directories ( 
  ${CMAKE_SOURCE_DIR}/include
  /usr/local/include 
  /opt/local/include 
  ${Boost_INCLUDE_DIRS} 
  ${USER_INCLUDE_PATH} 
  ${X11_X11_INCLUDE_PATH} )

#######################
# Library Directories #
#######################

link_directories ( 
  ${USER_LIBRARY_PATH} 
  /usr/local/lib 
  /opt/local/lib )

###########
# Outputs #
###########

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

###########
# Recurse #
###########

add_subdirectory (source)

###################
# Install Headers #
###################

file(GLOB HEADERS ${CMAKE_SOURCE_DIR}/include/chomp/*.h)
install(DIRECTORY include/ DESTINATION include/ )

#########
# Tests #
#########

enable_testing()
set(bin ${CMAKE_BINARY_DIR}/bin )
set(examples "${CMAKE_SOURCE_DIR}/examples")
add_test( circle ${CMAKE_SOURCE_DIR}/tests/circle.test ${bin} ${examples})
add_test( rubik ${CMAKE_SOURCE_DIR}/tests/segment.test ${bin} ${examples})
add_test( simplex ${CMAKE_SOURCE_DIR}/tests/simplex.test ${bin} ${examples})
add_test( square ${CMAKE_SOURCE_DIR}/tests/square.test ${bin} ${examples})
add_test( square_relative ${CMAKE_SOURCE_DIR}/tests/square_relative.test ${bin} ${examples})
add_test( square_withhole ${CMAKE_SOURCE_DIR}/tests/square_withhole.test ${bin} ${examples})

