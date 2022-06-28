file(GLOB HEADER_LIST CONFIGURE_DEPENDS "${MPS_SOURCE_DIR}/include/mps/*.h")
# set(HEADER_LIST "${MPS_SOURCE_DIR}/include/mps/*.h")

add_library(cuPC STATIC cuPC-S.cu ${HEADER_LIST})
target_include_directories(cuPC PUBLIC ../include)
target_compile_features(cuPC PUBLIC cxx_std_17)

add_library(corr STATIC corr_compressed.cu ${HEADER_LIST})
target_include_directories(corr PUBLIC ../include)
target_compile_features(corr PUBLIC cxx_std_17)

add_library(prep STATIC prep_markers.cpp ${HEADER_LIST})
set_source_files_properties(prep_markers.cpp PROPERTIES LANGUAGE CUDA)
target_include_directories(prep PUBLIC ../include)
target_compile_features(prep PUBLIC cxx_std_17)

add_library(io STATIC io.cpp ${HEADER_LIST})
target_include_directories(io PUBLIC ../include)
target_compile_features(io PUBLIC cxx_std_17)
