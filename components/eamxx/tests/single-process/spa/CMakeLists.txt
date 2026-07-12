include (ScreamUtils)

set (TEST_BASE_NAME spa_standalone)
set (FIXTURES_BASE_NAME ${TEST_BASE_NAME}_generate_output_nc_files)

# Create the test
CreateADUnitTest(${TEST_BASE_NAME}
  LABELS spa physics
  LIBS spa
  MPI_RANKS ${TEST_RANK_START} ${TEST_RANK_END}
  FIXTURES_SETUP_INDIVIDUAL ${FIXTURES_BASE_NAME}
)

# Set AD configurable options
SetVarDependingOnTestSize(NUM_STEPS 2 4 48)
set (ATM_TIME_STEP 864000)
set (RUN_T0 2021-10-12-45000)

## Copy (and configure) yaml files needed by tests
set (POSTFIX "")
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/input.yaml
               ${CMAKE_CURRENT_BINARY_DIR}/input.yaml)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/output.yaml
               ${CMAKE_CURRENT_BINARY_DIR}/output.yaml)

# Ensure test input files are present in the data dir
set (TEST_INPUT_FILES
  scream/init/${EAMxx_tests_IC_FILE_72lev}
  scream/maps/map_ne4np4_to_ne2np4_mono.nc
  scream/init/spa_file_unified_and_complete_ne4_20220428.nc
)
foreach (file IN ITEMS ${TEST_INPUT_FILES})
  GetInputFile(${file})
endforeach()

# Compare output files produced by npX tests, to ensure they are bfb
include (CompareNCFiles)

CompareNCFilesFamilyMpi (
  TEST_BASE_NAME ${TEST_BASE_NAME}
  FILE_META_NAME ${TEST_BASE_NAME}_output.INSTANT.nsteps_x1.npMPIRANKS.${RUN_T0}.nc
  MPI_RANKS ${TEST_RANK_START} ${TEST_RANK_END}
  LABELS spa physics
  META_FIXTURES_REQUIRED ${FIXTURES_BASE_NAME}_npMPIRANKS_omp1
)

if (SCREAM_ENABLE_BASELINE_TESTS)
  # Compare one of the output files with the baselines.
  # Note: one is enough, since we already check that np1 is BFB with npX
  set (OUT_FILE ${TEST_BASE_NAME}_output.INSTANT.nsteps_x1.np${TEST_RANK_END}.${RUN_T0}.nc)
  CreateBaselineTest(${TEST_BASE_NAME} ${TEST_RANK_END} ${OUT_FILE} ${FIXTURES_BASE_NAME})
endif()

if (EAMXX_ENABLE_PYTHON)
  #####################################
  #             SPA (JAX)             #
  #####################################
  # Swap the whole SPA step (time interpolation + vertical remap) for the
  # scream_jax implementation and compare against the C++ run.

  set (POSTFIX _jax)
  set (PY_MODULE_NAME "spa_jax")
  set (PY_MODULE_PATH ${SCREAM_BASE_DIR}/../../jax_port/scream_jax/adapters/eamxx)
  set (PY_BACKEND "host")
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/input.yaml
                 ${CMAKE_CURRENT_BINARY_DIR}/input_jax.yaml)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/output.yaml
                 ${CMAKE_CURRENT_BINARY_DIR}/output_jax.yaml)

  CreateUnitTestFromExec(${TEST_BASE_NAME}_jax ${TEST_BASE_NAME}
    EXE_ARGS "--args -ifile=input_jax.yaml"
    LABELS spa physics
    FIXTURES_SETUP spa_jax)

  # Tier-1 golden replay shows machine-precision agreement, so tight
  # tolerances apply.
  set (SRC_FILE "${TEST_BASE_NAME}_output.INSTANT.nsteps_x1.np1.${RUN_T0}.nc")
  set (TGT_FILE "${TEST_BASE_NAME}_jax_output.INSTANT.nsteps_x1.np1.${RUN_T0}.nc")
  add_test (NAME ${TEST_BASE_NAME}_cpp_vs_jax
            COMMAND ${Python_EXECUTABLE} ${SCREAM_BASE_DIR}/../../jax_port/harness/compare_nc_tolerant.py
            -s ${SRC_FILE} -t ${TGT_FILE} --tol 1e-10 --hard-tol 1e-8 --max-frac 1e-4
            -v nccn aero_g_sw aero_ssa_sw aero_tau_sw aero_tau_lw
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
  set_tests_properties(${TEST_BASE_NAME}_cpp_vs_jax PROPERTIES
        LABELS "spa;infrastructure"
        FIXTURES_REQUIRED "spa_jax;${FIXTURES_BASE_NAME}_np1_omp1")
endif()
