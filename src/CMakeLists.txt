file(GLOB HEADER_LIST CONFIGURE_DEPENDS "${MPS_SOURCE_DIR}/include/mps/*.h")
# set(HEADER_LIST "${MPS_SOURCE_DIR}/include/mps/*.h")

add_library(cuPC STATIC cuPC-S.cu cuPC_call_prep.cpp ${HEADER_LIST})
target_include_directories(cuPC PUBLIC ../include)
target_compile_features(cuPC PUBLIC cxx_std_20)

add_library(corr STATIC corr_host.cu corr_kernels.cu ${HEADER_LIST})
target_include_directories(corr PUBLIC ../include)
target_compile_features(corr PUBLIC cxx_std_20)

add_library(prep STATIC prep.cpp ${HEADER_LIST})
target_include_directories(prep PUBLIC ../include)
target_compile_features(prep PUBLIC cxx_std_20)

add_library(io STATIC io.cpp bim.cpp phen.cpp ${HEADER_LIST})
target_include_directories(io PUBLIC ../include)
target_compile_features(io PUBLIC cxx_std_20)

add_library(block STATIC blocking.cpp ${HEADER_LIST})
target_include_directories(block PUBLIC ../include)
target_compile_features(block PUBLIC cxx_std_20)

# add_library(cuPCsparse STATIC cu_test_helpers.cu cuPC-S-sparse.cu ${HEADER_LIST})
# target_include_directories(cuPCsparse PUBLIC ../include)
# target_compile_features(cuPCsparse PUBLIC cxx_std_20)


