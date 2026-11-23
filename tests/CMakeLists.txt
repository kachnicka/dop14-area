add_executable(dop14_tests tests.cpp)
target_link_libraries(dop14_tests PRIVATE dop14_core)

add_test(NAME dop14_tests COMMAND dop14_tests)
add_test(NAME jobs_match COMMAND "${CMAKE_COMMAND}" -DEXE=$<TARGET_FILE:dop14-area> -P "${CMAKE_CURRENT_SOURCE_DIR}/jobs_match.cmake")
