include (ScreamUtils)

set (TEST_BASE_NAME p3_standalone)
set (FIXTURES_BASE_NAME ${TEST_BASE_NAME}_generate_output_nc_files)

# Create the test
CreateADUnitTest(${TEST_BASE_NAME}
  LABELS p3 physics
  LIBS p3
  MPI_RANKS ${TEST_RANK_START} ${TEST_RANK_END}
  FIXTURES_SETUP_INDIVIDUAL ${FIXTURES_BASE_NAME}
)

# Set AD configurable options
set (ATM_TIME_STEP 1800)
SetVarDependingOnTestSize(NUM_STEPS 2 5 48)  # 1h 2.5h 24h
set (RUN_T0 2021-10-12-45000)

## Copy (and configure) yaml files needed by tests
set (POSTFIX "")
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/input.yaml
               ${CMAKE_CURRENT_BINARY_DIR}/input.yaml)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/output.yaml
               ${CMAKE_CURRENT_BINARY_DIR}/output.yaml)

# Ensure test input files are present in the data dir
GetInputFile(scream/init/${EAMxx_tests_IC_FILE_72lev})

# Compare output files produced by npX tests, to ensure they are bfb
include (CompareNCFiles)

CompareNCFilesFamilyMpi (
  TEST_BASE_NAME ${TEST_BASE_NAME}
  FILE_META_NAME ${TEST_BASE_NAME}_output.INSTANT.nsteps_x1.npMPIRANKS.${RUN_T0}.nc
  MPI_RANKS ${TEST_RANK_START} ${TEST_RANK_END}
  LABELS p3 physics
  META_FIXTURES_REQUIRED ${FIXTURES_BASE_NAME}_npMPIRANKS_omp1
)

# Check tendency calculation
foreach (NRANKS RANGE ${TEST_RANK_START} ${TEST_RANK_END})
  set (script ${SCREAM_BASE_DIR}/scripts/check-tendencies)
  set (fname p3_standalone_output.INSTANT.nsteps_x1.np${NRANKS}.${RUN_T0}.nc)
  set (tname p3_tend_check_np${NRANKS})
  add_test (NAME ${tname}
    COMMAND ${script} -f ${fname} -v qc T_mid -t p3_qc_tend p3_T_mid_tend
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
  set_tests_properties (${tname} PROPERTIES
    LABELS "p3;physics"
    FIXTURES_REQUIRED ${FIXTURES_BASE_NAME}_np${NRANKS}_omp1)
endforeach()

if (SCREAM_ENABLE_BASELINE_TESTS)
  # Compare one of the output files with the baselines.
  # Note: one is enough, since we already check that np1 is BFB with npX
  set (OUT_FILE ${TEST_BASE_NAME}_output.INSTANT.nsteps_x1.np${TEST_RANK_END}.${RUN_T0}.nc)
  CreateBaselineTest(${TEST_BASE_NAME} ${TEST_RANK_END} ${OUT_FILE} ${FIXTURES_BASE_NAME}) 

  # Grab all bfb hash lines into a text file (only np1_omp1 test)
  CreateUnitTest(${TEST_BASE_NAME}_bfb_hash
    SOURCES p3_bfb_hash.cpp
    FIXTURES_REQUIRED ${FIXTURES_BASE_NAME}_np1_omp1
    FIXTURES_SETUP ${TEST_BASE_NAME}_gen_bfb_hash
    EXE_ARGS "--args --log-file atm.log.np1 --output-file p3_bfb_hash.txt"
    LABELS p3 baseline_gen
  )

  # Compare against bfb hashes in baseline dir
  set (SRC_FILE ${CMAKE_CURRENT_BINARY_DIR}/p3_bfb_hash.txt)
  set (TGT_FILE ${SCREAM_BASELINES_DIR}/data/p3_bfb_hash.txt)
  add_test (NAME ${TEST_BASE_NAME}_bfb_hash_baseline_cmp
    COMMAND ${CMAKE_COMMAND} -P ${SCREAM_BASE_DIR}/cmake/CompareTextFiles.cmake ${SRC_FILE} ${TGT_FILE})

  set_tests_properties (
    ${TEST_BASE_NAME}_bfb_hash_baseline_cmp PROPERTIES
    FIXTURES_REQUIRED ${TEST_BASE_NAME}_gen_bfb_hash
    LABELS "p3;baseline_cmp")

  # Store p3_bfb_hash.txt in the baseline folder
  file (APPEND ${SCREAM_TEST_OUTPUT_DIR}/baseline_list
    "${SRC_FILE}\n"
  )
endif()

if (EAMXX_ENABLE_PYTHON)
  #####################################
  #             P3 (JAX)              #
  #####################################
  # Swap the whole P3 step for the scream_jax implementation (see
  # jax_port/TEST_HARNESS_DESIGN.md) and compare against the C++ run.
  # Requires `jax` importable by the embedded interpreter.

  set (POSTFIX _jax)
  set (PY_MODULE_NAME "p3_jax")
  set (PY_MODULE_PATH ${SCREAM_BASE_DIR}/../../jax_port/scream_jax/adapters/eamxx)
  set (PY_BACKEND "host")
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/input.yaml
                 ${CMAKE_CURRENT_BINARY_DIR}/input_jax.yaml)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/output.yaml
                 ${CMAKE_CURRENT_BINARY_DIR}/output_jax.yaml)

  CreateUnitTestFromExec(${TEST_BASE_NAME}_jax ${TEST_BASE_NAME}
    EXE_ARGS "--args -ifile=input_jax.yaml"
    LABELS p3 physics
    FIXTURES_SETUP p3_jax)

  # Compare against the C++ (np1) run. Tolerance-based (the JAX port is
  # not BFB). Single-step agreement is at roundoff level for the whole
  # warm/thermo path (see jax_port/tests/test_p3_golden.py), but ice
  # sedimentation amplifies roundoff through discontinuous branches
  # (verified inherent by feeding the REAL C++ ice sed two part2 states
  # differing by 1e-13: same divergence), and over this multi-step
  # coupled run that noise feeds back into the liquid fields (melt/
  # collection). T_mid/qv stay roundoff-clean over all steps and keep
  # tight bounds; everything touching the ice feedback gets a
  # bounded-outlier criterion.
  set (SRC_FILE "${TEST_BASE_NAME}_output.INSTANT.nsteps_x1.np1.${RUN_T0}.nc")
  set (TGT_FILE "${TEST_BASE_NAME}_jax_output.INSTANT.nsteps_x1.np1.${RUN_T0}.nc")
  add_test (NAME ${TEST_BASE_NAME}_cpp_vs_jax_strict
            COMMAND ${Python_EXECUTABLE} ${SCREAM_BASE_DIR}/../../jax_port/harness/compare_nc_tolerant.py
            -s ${SRC_FILE} -t ${TGT_FILE} --tol 1e-6 --hard-tol 1e-2 --max-frac 1e-3
            -v T_mid T_prev_micro_step qv qv_prev_micro_step rainfrac
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
  set_tests_properties(${TEST_BASE_NAME}_cpp_vs_jax_strict PROPERTIES
        LABELS "p3;infrastructure"
        FIXTURES_REQUIRED "p3_jax;${FIXTURES_BASE_NAME}_np1_omp1")
  add_test (NAME ${TEST_BASE_NAME}_cpp_vs_jax_ice_feedback
            COMMAND ${Python_EXECUTABLE} ${SCREAM_BASE_DIR}/../../jax_port/harness/compare_nc_tolerant.py
            -s ${SRC_FILE} -t ${TGT_FILE} --tol 1e-6 --hard-tol 0.5 --max-frac 0.1
            -v qc nc qr nr qi ni qm bm
               eff_radius_qc eff_radius_qr eff_radius_qi
               precip_liq_surf_mass precip_ice_surf_mass
               micro_liq_ice_exchange micro_vap_liq_exchange
               micro_vap_ice_exchange
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
  set_tests_properties(${TEST_BASE_NAME}_cpp_vs_jax_ice_feedback PROPERTIES
        LABELS "p3;infrastructure"
        FIXTURES_REQUIRED "p3_jax;${FIXTURES_BASE_NAME}_np1_omp1")
endif()
