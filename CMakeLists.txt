cmake_minimum_required(VERSION 3.2)

include(FetchContent)
FetchContent_Declare(
  vcpkg
  GIT_REPOSITORY https://github.com/microsoft/vcpkg
)
FetchContent_MakeAvailable(vcpkg)
FetchContent_GetProperties(vcpkg)

message(STATUS ${vcpkg_SOURCE_DIR})
if(NOT CMAKE_TOOLCHAIN_FILE)
    set(CMAKE_TOOLCHAIN_FILE "${vcpkg_SOURCE_DIR}\\scripts\\buildsystems\\vcpkg.cmake")
endif()
if(NOT VCPKG_MANIFEST_MODE)
    set(VCPKG_MANIFEST_MODE ON)
endif()

project(dsp)


find_package(FFTW3 REQUIRED)
find_package(Matplot++ REQUIRED)
set(SRCS    "src/cx_math.cpp"
            "src/naive.cpp"
            "src/viz.cpp"
            "src/utils.cpp"
            "src/main.cpp")

add_executable(dsp ${SRCS})
target_include_directories(dsp PUBLIC "include")

SET(LIBS FFTW3::fftw3 Matplot++::matplot)
target_link_libraries(dsp PUBLIC ${LIBS})

foreach(LIB ${LIBS})
    get_target_property(INC_DIR ${LIB} INTERFACE_INCLUDE_DIRECTORIES)
    target_include_directories(dsp PUBLIC ${INC_DIR})
endforeach()