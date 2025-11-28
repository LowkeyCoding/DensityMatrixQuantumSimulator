# Add source files
add_library(dmqs_core SHARED
    dmqs.cpp
    channels.cpp
    gates.cpp
)

target_include_directories(dmqs_core PUBLIC 
    ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(dmqs_core PUBLIC 
    armadillo 
    openblas)