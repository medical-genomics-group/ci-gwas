file(GLOB HEADER_LIST CONFIGURE_DEPENDS "${MPS_SOURCE_DIR}/include/mps/*.h")

# set(HEADER_LIST "${MPS_SOURCE_DIR}/include/mps/*.h")
add_library(cuPC STATIC cuPC-S.cu hetcor-cuPC-S.cu cuPC-S-second-stage.cu cuPC_call_prep.cpp ${HEADER_LIST})
set_target_properties(cuPC PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
target_include_directories(cuPC PUBLIC ../include)
target_compile_features(cuPC PUBLIC cxx_std_20)
target_compile_options(cuPC PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:
    -dc
    -O3
    -lineinfo
    --use_fast_math
    --shared
    -Xcompiler
    -fPIC
    >)

add_library(corr STATIC corr_host.cu corr_kernels.cu ${HEADER_LIST})
target_include_directories(corr PUBLIC ../include)
target_compile_features(corr PUBLIC cxx_std_20)
target_compile_options(corr PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:
    -O3
    -lineinfo
    --use_fast_math
    --shared
    -Xcompiler
    -fPIC
    >)

add_library(prep STATIC prep.cpp ${HEADER_LIST})
target_include_directories(prep PUBLIC ../include)
target_compile_features(prep PUBLIC cxx_std_20)

add_library(io STATIC io.cpp bim.cpp phen.cpp trait_summary_stats.cpp marker_trait_summary_stats.cpp marker_summary_stats.cpp ${HEADER_LIST})
target_include_directories(io PUBLIC ../include)
target_compile_features(io PUBLIC cxx_std_20)

add_library(block STATIC blocking.cpp ${HEADER_LIST})
target_include_directories(block PUBLIC ../include)
target_compile_features(block PUBLIC cxx_std_20)

add_library(parents STATIC parent_set.cpp ${HEADER_LIST})
target_include_directories(parents PUBLIC ../include)
target_compile_features(parents PUBLIC cxx_std_20)
