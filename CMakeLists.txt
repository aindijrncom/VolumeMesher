cmake_minimum_required(VERSION 3.10)

# Set the target architecture.
set(ENABLE_SSE2 True)
set(ENABLE_AVX2 True)

#####################

# use gcc and g++ instead of clang (native build only)
if(NOT CMAKE_CROSSCOMPILING)
    set(CMAKE_CXX_COMPILER "g++-12")
    set(CMAKE_C_COMPILER "gcc-12")
endif()

# set the project name
project(mesh_generator)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set (SOURCES
    src/makePolyhedralMesh.cpp
    src/delaunay.cpp
    src/conforming_mesh.cpp
    src/extended_predicates.cpp
    src/BSP.cpp
    src/inOutPartition.cpp
    src/winMem.cpp
)

# add the executable
add_executable(${PROJECT_NAME}
    src/main.cpp
    ${SOURCES}
)

# make a static library for use within external programs
add_library(${PROJECT_NAME}_lib
    ${SOURCES}
)

set(ALL_TARGETS ${PROJECT_NAME} ${PROJECT_NAME}_lib)

foreach (TARGET ${ALL_TARGETS})
if(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
	# grant IEEE 754 compliance
	target_compile_options(${TARGET} PUBLIC "/fp:strict")
	# use intrinsic functions
	target_compile_options(${TARGET} PUBLIC "/Oi")
	# enable AVX2 instruction set
	target_compile_options(${TARGET} PUBLIC "/arch:AVX2")
	# turn off annoying warnings
	target_compile_options(${TARGET} PUBLIC "/D _CRT_SECURE_NO_WARNINGS")
	# set target architecture
	if(ENABLE_AVX2)
		target_compile_options(${PROJECT_NAME} PUBLIC "/arch:AVX2")
	elseif(ENABLE_SSE2)
		target_compile_options(${PROJECT_NAME} PUBLIC "/arch:SSE2")
	endif()
else()
	# set standard optimization level
	target_compile_options(${TARGET} PUBLIC -O2)
	# grant IEEE 754 compliance
	target_compile_options(${TARGET} PUBLIC -frounding-math)

	# set target architecture
	if(ENABLE_AVX2)
		target_compile_options(${PROJECT_NAME} PUBLIC "-mavx2")
		target_compile_options(${PROJECT_NAME} PUBLIC "-mfma")
	elseif(ENABLE_SSE2)
		target_compile_options(${PROJECT_NAME} PUBLIC "-msse2")
	endif()

	# ---- Release hardening: reduce information leaks ----
	# Hide all symbols (export table reveals class/function names)
	target_compile_options(${TARGET} PRIVATE -fvisibility=hidden -fvisibility-inlines-hidden)
	# Remove GCC version string from .comment section
	target_compile_options(${TARGET} PRIVATE -fno-ident)
	# Reduce .eh_frame path leaks (fewer source paths in binary)
	# Skip on MinGW: interferes with __gxx_personality_v0 (C++ exception handling)
	if(NOT CMAKE_SYSTEM_NAME STREQUAL "Windows")
		target_compile_options(${TARGET} PRIVATE $<$<CONFIG:Release>:-fno-asynchronous-unwind-tables>)
	endif()

	# Strip executable in Release mode
	if(CMAKE_BUILD_TYPE STREQUAL "Release")
		target_link_options(${TARGET} PRIVATE -Wl,-s)
	endif()

	# MinGW cross-compilation: static link to minimize DLL dependencies
	if(CMAKE_SYSTEM_NAME STREQUAL "Windows")
		target_link_options(${TARGET} PRIVATE
			-static -static-libgcc -static-libstdc++
			-Wl,--allow-multiple-definition)
	endif()
endif()

# Public include directory
target_include_directories(${TARGET} PUBLIC
	src
)
endforeach(TARGET)
