

if(GATES_TESTS)
    add_executable(gates_test gates_test.cpp)
    target_link_libraries(gates_test dmqs_core doctest::doctest_with_main)

    add_custom_command(
        TARGET gates_test
        POST_BUILD
        COMMAND $<TARGET_FILE:gates_test>
    )
endif(GATES_TESTS)

if(DMQS_TESTS)
    add_executable(dmqs_test dmqs_test.cpp)
    target_link_libraries(dmqs_test dmqs_core doctest::doctest_with_main)

    add_custom_command(
        TARGET dmqs_test
        POST_BUILD
        COMMAND $<TARGET_FILE:dmqs_test>
    )
endif(DMQS_TESTS)

if(UPPAAL_TESTS)
    add_executable(uppaal_test uppaal_test.cpp)
    target_link_libraries(uppaal_test dmqs_uppaal doctest::doctest_with_main)

    add_custom_command(
        TARGET uppaal_test
        POST_BUILD
        COMMAND $<TARGET_FILE:uppaal_test>
    )
endif(UPPAAL_TESTS)