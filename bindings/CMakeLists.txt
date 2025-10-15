# Create the shared library with uppaal binding
add_library(dmqs_uppaal SHARED
    uppaal.cpp
)

target_link_libraries(dmqs_uppaal PUBLIC 
    dmqs_core
)