cmake_minimum_required(VERSION 3.10)
if(POLICY CMP0135)
    cmake_policy(SET CMP0135 NEW)
endif()

# Set the target architecture.
set(ENABLE_SSE2 True)
set(ENABLE_AVX2 True)
# set(ENABLE_ARMNEON True)

if(APPLE)
    set(ENABLE_SSE2 False)
    set(ENABLE_AVX2 False)
endif()

if(ENABLE_AVX2)
    set(ENABLE_SSE2 True)
endif()

if(ENABLE_ARMNEON AND ENABLE_SSE2)
    message(FATAL_ERROR "ARMNEON and SSE2/AVX2 cannot be both set!")
endif()

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

# Upstream geometric-number dependencies.
include(FetchContent)

if(EXISTS ../.local.cmake)
    set(nfg_SOURCE_DIR "../NFG/")
    set(indirect_predicates_SOURCE_DIR "../Indirect_Predicates/")
else()
    FetchContent_Declare(
        nfg
        URL https://github.com/MarcoAttene/nfg/archive/refs/heads/master.zip)
    FetchContent_Populate(nfg)
    message(STATUS "NFG is available in " ${nfg_SOURCE_DIR})

    FetchContent_Declare(
        Indirect_Predicates
        URL https://github.com/MarcoAttene/Indirect_Predicates/archive/refs/heads/master.zip)
    FetchContent_Populate(Indirect_Predicates)
    message(STATUS "Indirect_Predicates is available in " ${indirect_predicates_SOURCE_DIR})
endif()

if(ENABLE_ARMNEON)
    FetchContent_Declare(
        simde
        URL https://github.com/simd-everywhere/simde/archive/refs/heads/master.zip)
    FetchContent_Populate(simde)
    message(STATUS "simde is available in " ${simde_SOURCE_DIR})
    target_include_directories(${PROJECT_NAME} PUBLIC ${simde_SOURCE_DIR}/simde)
    target_compile_options(${PROJECT_NAME} PUBLIC -DSIMDE_ENABLE_NATIVE_ALIASES)
endif()

set(ALL_TARGETS ${PROJECT_NAME} ${PROJECT_NAME}_lib)

foreach (TARGET ${ALL_TARGETS})
if(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
	# grant IEEE 754 compliance
	target_compile_options(${TARGET} PUBLIC "/fp:strict")
	# use intrinsic functions
	target_compile_options(${TARGET} PUBLIC "/Oi")
	# turn off annoying warnings
	target_compile_options(${TARGET} PUBLIC "/D _CRT_SECURE_NO_WARNINGS")
	# set target architecture
	if(ENABLE_AVX2)
		target_compile_options(${TARGET} PUBLIC "/arch:AVX2")
	elseif(ENABLE_SSE2)
		target_compile_options(${TARGET} PUBLIC "/arch:SSE2")
	elseif(ENABLE_ARMNEON)
		target_compile_options(${TARGET} PUBLIC "/arch:ARMV8.0")
	endif()
else()
	# set standard optimization level
	target_compile_options(${TARGET} PUBLIC -O3)
	# grant IEEE 754 compliance
	target_compile_options(${TARGET} PUBLIC -frounding-math)
	# Preserve exact-predicate error-free transformations.
	target_compile_options(${TARGET} PUBLIC -ffp-contract=off)

	# set target architecture
	if(ENABLE_AVX2)
		target_compile_options(${TARGET} PUBLIC "-mavx2")
		target_compile_options(${TARGET} PUBLIC "-mfma")
	elseif(ENABLE_SSE2)
		target_compile_options(${TARGET} PUBLIC "-msse2")
	elseif(ENABLE_ARMNEON)
		target_compile_options(${TARGET} PUBLIC "-march=armv8-a+simd")
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
	${nfg_SOURCE_DIR}/include
	${indirect_predicates_SOURCE_DIR}/include
)
endforeach(TARGET)
