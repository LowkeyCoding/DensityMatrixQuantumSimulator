

if(GATES_TESTS)
    add_executable(gates_test gates_test.cpp)
    target_link_libraries(gates_test dmqs_core)

    # Add test running command
    add_custom_command(
        TARGET gates_test
        POST_BUILD
        COMMAND $<TARGET_FILE:gates_test>
    )
endif(GATES_TESTS)

if(DMQS_TESTS)
    add_executable(dmqs_test dmqs_test.cpp)
    target_link_libraries(dmqs_test dmqs_core)

    # Add test running command
    add_custom_command(
        TARGET dmqs_test
        POST_BUILD
        COMMAND $<TARGET_FILE:dmqs_test>
    )
endif(DMQS_TESTS)