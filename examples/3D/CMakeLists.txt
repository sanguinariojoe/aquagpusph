ADD_SUBDIRECTORY(spheric_testcase2_dambreak)
ADD_SUBDIRECTORY(spheric_testcase9_tld)
ADD_SUBDIRECTORY(spheric_testcase10_waveimpact)

IF(HAVE_MPI)
    ADD_SUBDIRECTORY(spheric_testcase2_dambreak_mpi)
ENDIF(HAVE_MPI)
