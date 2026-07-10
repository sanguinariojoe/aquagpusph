add_subdirectory(spheric_testcase2_dambreak)
add_subdirectory(spheric_testcase9_tld)
add_subdirectory(spheric_testcase10_waveimpact)

if(HAVE_MPI)
    add_subdirectory(spheric_testcase2_dambreak_mpi)
endif(HAVE_MPI)

if(HAVE_CHRONO)
    add_subdirectory(apollo_capsule)
    add_subdirectory(flow_through_rocks)
    add_subdirectory(shock_rocks_BIe)
    add_subdirectory(ball_inflow_britan_etal_1995_midpoint)
endif(HAVE_CHRONO)
