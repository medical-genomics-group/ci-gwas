add_executable(main main.cpp)
set_source_files_properties(main.cpp PROPERTIES LANGUAGE CUDA)
target_compile_features(main PRIVATE cxx_std_17)

target_link_libraries(main PRIVATE corr)
