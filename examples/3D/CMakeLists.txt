add_subdirectory(spheric_testcase2_dambreak)
add_subdirectory(spheric_testcase9_tld)
add_subdirectory(spheric_testcase10_waveimpact)

if(HAVE_MPI)
    add_subdirectory(spheric_testcase2_dambreak_mpi)
endif(HAVE_MPI)

if(HAVE_CHRONO)
    add_subdirectory(titan_capsule)
endif(HAVE_CHRONO)
