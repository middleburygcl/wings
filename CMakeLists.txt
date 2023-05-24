cmake_minimum_required(VERSION 3.1)

project(wings LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "-O3")

set(LIBRARIES)
if (APPLE)
	find_library(OpenGL_LIBRARY OpenGL)
	set(LIBRARIES ${LIBRARIES} ${OpenGL_LIBRARY})
else()
	find_package(OpenGL COMPONENTS REQUIRED OpenGL EGL)
	set(LIBRARIES ${LIBRARIES} OpenGL::EGL pthread)
endif()

function(add_extern_repository name)
	set(options "")
	set(one_value_args GIT_REPOSITORY FULL_HISTORY SKIP_CONFIG)
	set(multi_value_args "")
	cmake_parse_arguments(ARG "${options}" "${one_value_args}" "${multi_value_args}" ${ARGN})
	if(NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/extern/${name})	
		if(ARG_GIT_REPOSITORY)
			set(clone_opts --recursive --depth=1 --single-branch)
			if(ARG_FULL_HISTORY)
				set(clone_opts --recursive)
			endif()
			set(fetch_get git)
			set(fetch_url ${ARG_GIT_REPOSITORY})
			set(fetch_arg clone ${clone_opts} ${ARG_GIT_REPOSITORY} ${CMAKE_CURRENT_SOURCE_DIR}/${name})
		else()
			message(FATAL_ERROR "unknown repository type")
		endif()
		message(STATUS "fetching ${name} from ${fetch_url}")
		execute_process(COMMAND ${fetch_get} ${fetch_arg} RESULT_VARIABLE status OUTPUT_QUIET ERROR_QUIET)
	endif()
endfunction()

add_library(wings SHARED wings.cpp)
set_property(TARGET wings PROPERTY CXX_STANDARD 14)
target_link_libraries(wings ${LIBRARIES})

option(BUILD_EXAMPLE "build example" TRUE)
set(WINGS_STB_LOCATION "stb" CACHE STRING "location of stb installation")
if (${WINGS_STB_LOCATION} MATCHES "stb")
	add_extern_repository(stb GIT_REPOSITORY "https://github.com/nothings/stb")
	set(WINGS_STB_LOCATION ${CMAKE_CURRENT_SOURCE_DIR}/stb)
	set_directory_properties(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES ${CMAKE_CURRENT_SOURCE_DIR}/stb)
	add_definitions(-DWINGS_COMPILE_STB=1)
endif()

include_directories(${WINGS_STB_LOCATION})
if (${BUILD_EXAMPLE})
	add_extern_repository(tinyobjloader GIT_REPOSITORY "https://github.com/tinyobjloader/tinyobjloader")
	add_executable(wings_example example/main.cpp tinyobjloader/tiny_obj_loader.cc)
	target_link_libraries(wings_example wings)
	set_property(TARGET wings_example PROPERTY CXX_STANDARD 14)
endif()

set_property(DIRECTORY APPEND PROPERTY ADDITIONAL_MAKE_CLEAN_FILES 
	${CMAKE_CURRENT_SOURCE_DIR}/stb
	${CMAKE_CURRENT_SOURCE_DIR}/tinyobjloader
)



