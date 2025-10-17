# Ensures that Armadillo (Linear Algebra) library is installed
include(FetchContent)
FetchContent_Declare(Armadillo
        URL https://sourceforge.net/projects/arma/files/armadillo-15.0.3.tar.xz
        URL_HASH SHA256=9f55ec10f0a91fb6479ab4ed2b37a52445aee917706a238d170b5220c022fe43
        DOWNLOAD_NO_PROGRESS ON
        EXCLUDE_FROM_ALL ON
        FIND_PACKAGE_ARGS 15)
set(STATIC_LIB ON CACHE BOOL "Build static Armadillo library")
set(OPENBLAS_PROVIDES_LAPACK ON CACHE BOOL "Assume that OpenBLAS provides LAPACK functions")
set(ALLOW_FLEXIBLAS_LINUX ON CACHE BOOL "Allow detection of FlexiBLAS on Linux")
set(ALLOW_OPENBLAS_MACOS ON CACHE BOOL "Allow detection of OpenBLAS on macOS")
set(ALLOW_BLAS_LAPACK_MACOS ON CACHE BOOL "Allow detection of BLAS and LAPACK on macOS")
set(BUILD_SMOKE_TEST OFF CACHE BOOL "Build the Armadillo smoke test")
FetchContent_MakeAvailable(Armadillo)
if (Armadillo_FOUND) # by find_package
    message(STATUS "Found Armadillo: ${ARMADILLO_INCLUDE_DIRS}")
else (Armadillo_FOUND) # by FetchContent
    message(STATUS "Fetched Armadillo: ${armadillo_SOURCE_DIR}")
endif (Armadillo_FOUND)
if (TARGET Armadillo::Armadillo)
    message(STATUS "    Available target: Armadillo::Armadillo")
else ()
    if (TARGET armadillo)
        message(STATUS "    Available target: armadillo")
        add_library(Armadillo::Armadillo INTERFACE IMPORTED)
        set_target_properties(Armadillo::Armadillo
                PROPERTIES
                INTERFACE_LINK_LIBRARIES armadillo)
        message(STATUS "    Added target: Armadillo::Armadillo")
    else()
        message(STATUS "    ARMADILLO_INCLUDE_DIRS=${ARMADILLO_INCLUDE_DIRS}")
        message(STATUS "    ARMADILLO_LIBRARY_DIRS=${ARMADILLO_LIBRARY_DIRS}")
        message(STATUS "    ARMADILLO_LIBRARIES=${ARMADILLO_LIBRARIES}")
        add_library(Armadillo::Armadillo INTERFACE IMPORTED)
        set_target_properties(Armadillo::Armadillo
                PROPERTIES
                INTERFACE_LINK_LIBRARIES "${ARMADILLO_LIBRARIES}"
                INTERFACE_LINK_DIRECTORIES "${ARMADILLO_LIBRARY_DIRS}"
                INTERFACE_INCLUDE_DIRECTORIES "${ARMADILLO_INCLUDE_DIRS}")
        message(STATUS "    Added target: Armadillo::Armadillo")
    endif()
endif()
