
IF (NOT QMC_CUDA)

  #
  # H2 test in a DZV basis set, using optimizable determinants, and the adaptive linear method
  #
  LIST(APPEND HE_PARAM_SCALARS "jud_0" -0.122201  0.0065) # bspline param gradient
  LIST(APPEND HE_PARAM_SCALARS "jud_1"  0.06973   0.0077)
  LIST(APPEND HE_PARAM_SCALARS "jud_2"  0.05146   0.0044)
  LIST(APPEND HE_PARAM_SCALARS "jud_3"  0.000814  0.00036)

  MESSAGE_VERBOSE("Adding He-param tests")

  QMC_RUN_AND_CHECK2(short-He-param
                    "${CMAKE_SOURCE_DIR}/tests/molecules/He_param"
                    He_param_grad.param
                    He_param_grad.xml
                    1 1
                    TRUE 
                    ${HE_PARAM_SCALARS}
                    )

ELSE()
    MESSAGE_VERBOSE("Skipping He tests because gaussian basis sets are not supported by CUDA build (QMC_CUDA=1)")
ENDIF()

