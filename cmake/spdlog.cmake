# Ensures that Armadillo (Linear Algebra) library is installed
include(FetchContent)
FetchContent_Declare(spdlog
        GIT_REPOSITORY https://github.com/gabime/spdlog.git
        GIT_TAG v1.15.3   # "main" for latest
        GIT_SHALLOW TRUE  # download specific revision only (git clone --depth 1)
        GIT_PROGRESS TRUE # show download progress in Ninja
        EXCLUDE_FROM_ALL ON # don't build if not used
        FIND_PACKAGE_ARGS 1.15)
set(SPDLOG_BUILD_ALL OFF CACHE BOOL "Build all artifacts")
set(SPDLOG_BUILD_SHARED OFF CACHE BOOL "Build shared library")
set(SPDLOG_ENABLE_PCH OFF CACHE BOOL "Build precompiled headers to speedup compilation")
set(SPDLOG_BUILD_EXAMPLE OFF CACHE BOOL "Build example")
set(SPDLOG_BUILD_EXAMPLE_HO OFF CACHE BOOL "Build header-only example")
set(SPDLOG_BUILD_TESTS OFF CACHE BOOL "Build tests")
set(SPDLOG_BUILD_TESTS_HO OFF CACHE BOOL "Build header-only tests")
set(SPDLOG_BUILD_BENCH OFF CACHE BOOL "Build benchmarks (assumes GoogleBenchmark)")
set(SPDLOG_SANITIZE_ADDRESS OFF CACHE BOOL "Enable AddressSanitizer in tests")
set(SPDLOG_SANITIZE_THREAD OFF CACHE BOOL "Enable ThreadSanitizer in tests")
set(SPDLOG_BUILD_WARNINGS OFF CACHE BOOL "Enable compiler warnings")
set(SPDLOG_SYSTEM_INCLUDES OFF CACHE BOOL "Include as system headers (skip for clang-tidy).")
set(SPDLOG_INSTALL OFF CACHE BOOL "Generate install targets")
set(SPDLOG_USE_STD_FORMAT ON CACHE BOOL "Use std::format instead of fmt")
set(SPDLOG_FMT_EXTERNAL OFF CACHE BOOL "Use external fmd library instead of bundled")
set(SPDLOG_FMT_EXTERNAL_HO OFF CACHE BOOL "Use external fmd header-only library instead of bundled")
set(SPDLOG_NO_EXCEPTIONS OFF CACHE BOOL "Compile with -fno-exceptions. Call abort() on spdlog exceptions")
set(SPDLOG_PREVENT_CHILD_FD OFF CACHE BOOL "Prevent child processes from inheriting log file descriptors")
set(SPDLOG_NO_EXCEPTIONS OFF CACHE BOOL "prevent spdlog from querying the thread id on each log call if thread id is not needed")
set(SPDLOG_NO_TLS OFF CACHE BOOL "Prevent spdlog from using thread_local storage")
set(SPDLOG_NO_ATOMIC_LEVELS OFF CACHE BOOL "Prevent spdlog from using std::atomic log levels (only if your code does not modify log levels concurrently)")
set(SPDLOG_DISABLE_DEFAULT_LOGGER OFF CACHE BOOL "Disable creation of the default logger")
set(SPDLOG_FWRITE_UNLOCKED ON CACHE BOOL "Use the unlocked variant of fwrite. Leave this on unless your libc doesn't have it")
set(SPDLOG_TIDY OFF CACHE BOOL "Run clang-tody")

FetchContent_MakeAvailable(spdlog)
if (spdlog_FOUND) # by find_package
    message(STATUS "Found spdlog: ${spdlog_INCLUDE_DIRS}")
else (spdlog_FOUND) # by FetchContent
    message(STATUS "Fetched spdlog: ${spdlog_SOURCE_DIR}")
endif (spdlog_FOUND)
if (TARGET spdlog::spdlog)
    message(STATUS "    Available target: spdlog::spdlog")
endif()
