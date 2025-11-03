#!/bin/sh
current_branch=$(git symbolic-ref --short HEAD 2>/dev/null)
RUN_HOOK=true

if [[ "$RUN_HOOK" = true ]]; then 
    cmake -B build-release -DCMAKE_BUILD_TYPE=Release
    cmake --build build-release -- -j 8
    ctest --test-dir build-release --output-on-failure
else 
    echo "Skipping unit tests."
    exit 1
fi