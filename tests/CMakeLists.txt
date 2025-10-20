add_executable(gates_test gates_test.cpp)
target_link_libraries(gates_test dmqs_core doctest::doctest_with_main)
add_test(gates_test gates_test)


add_executable(dmqs_test dmqs_test.cpp)
target_link_libraries(dmqs_test dmqs_core doctest::doctest_with_main)
add_test(dmqs_test dmqs_test)

add_executable(uppaal_test uppaal_test.cpp)
target_link_libraries(uppaal_test dmqs_uppaal doctest::doctest_with_main)
add_test(uppaal_test uppaal_test)
