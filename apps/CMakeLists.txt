# add_executable(main main.cpp)
# set_source_files_properties(main.cpp PROPERTIES LANGUAGE CUDA)
# target_compile_features(main PRIVATE cxx_std_20)
# target_link_libraries(main PRIVATE corr)

add_executable(mps mps.cpp)
set_source_files_properties(mps.cpp PROPERTIES LANGUAGE CUDA)
target_compile_features(mps PRIVATE cxx_std_20)

# target_link_libraries(mps PRIVATE corr prep io cuPCsparse)
target_link_libraries(mps PRIVATE corr prep io cuPC parents block)
target_compile_options(mps PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:
    --generate-line-info
    --use_fast_math
    --relocatable-device-code=true
    >)