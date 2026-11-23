# cmake -DEXE=<dop14-area> -P jobs_match.cmake: a forked run's correctness
# report must match one process's, apart from its elapsed time; the timing
# table after it is not compared
foreach (jobs 1 3)
    execute_process(COMMAND "${EXE}" --count 150 --seed 11 --jobs ${jobs}
        OUTPUT_VARIABLE out RESULT_VARIABLE rc)
    if (NOT rc EQUAL 0)
        message(FATAL_ERROR "--jobs ${jobs}: exit code ${rc}")
    endif()
    string(FIND "${out}" "\n" begin)
    string(FIND "${out}" "time per DOP" end)
    math(EXPR length "${end} - ${begin}")
    string(SUBSTRING "${out}" ${begin} ${length} report_${jobs})
endforeach()
if (NOT report_1 STREQUAL report_3)
    message(FATAL_ERROR "--jobs 3 differs from --jobs 1")
endif()
