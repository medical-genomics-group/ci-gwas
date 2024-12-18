file(GLOB HEADER_LIST CONFIGURE_DEPENDS "${MPS_SOURCE_DIR}/include/mps/*.h")

add_library(mps_headers INTERFACE)
target_include_directories(mps_headers INTERFACE 
    ../include
    ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES}  # Important: Add CUDA include path
)
target_sources(mps_headers INTERFACE ${HEADER_LIST})

# # Force CUDA compilation for files with CUDA headers
# set_source_files_properties(
#     cli.cpp io.cpp prep.cpp block.cpp parents.cpp
#     PROPERTIES LANGUAGE CUDA
# )

add_library(cuPC STATIC
    cuPC-S.cu
    hetcor-cuPC-S.cu
    cuPC-S-second-stage.cu
    cuPC_call_prep.cpp
)
set_target_properties(cuPC PROPERTIES 
    CUDA_SEPARABLE_COMPILATION ON
    LINKER_LANGUAGE CUDA
)
target_link_libraries(cuPC PUBLIC mps_headers)
target_compile_features(cuPC PUBLIC cxx_std_20)
target_compile_options(cuPC PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:
    -dc
    -O3
    --expt-relaxed-constexpr
    -lineinfo
    --use_fast_math
    --shared
    -Xcompiler
    -fPIC
    >)

add_library(corr STATIC
    corr_host.cu
    corr_kernels.cu
)
target_link_libraries(corr PUBLIC mps_headers)
set_target_properties(corr PROPERTIES 
    CUDA_SEPARABLE_COMPILATION ON
    LINKER_LANGUAGE CUDA
)
target_compile_features(corr PUBLIC cxx_std_20)
target_compile_options(corr PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:
    -O3
    -lineinfo
    --use_fast_math
    --shared
    -Xcompiler
    -fPIC
    >)

add_library(prep STATIC
    prep.cpp
)
target_link_libraries(prep PUBLIC mps_headers)
target_compile_features(prep PUBLIC cxx_std_20)

add_library(io STATIC
    io.cpp
    bim.cpp
    phen.cpp
    trait_summary_stats.cpp
    marker_trait_summary_stats.cpp
    marker_summary_stats.cpp
)
target_link_libraries(io PUBLIC mps_headers)
target_compile_features(io PUBLIC cxx_std_20)

add_library(block STATIC
    blocking.cpp
)
target_link_libraries(block PUBLIC mps_headers)
target_compile_features(block PUBLIC cxx_std_20)

add_library(parents STATIC
    parent_set.cpp
)
target_link_libraries(parents PUBLIC mps_headers)
target_compile_features(parents PUBLIC cxx_std_20)

add_library(cli STATIC
    cli.cpp
)
set_target_properties(cli PROPERTIES 
    LINKER_LANGUAGE CUDA  # This is key
)
target_link_libraries(cli PUBLIC mps_headers io corr)
target_compile_features(cli PUBLIC cxx_std_20)