#pragma once
#include <string>
#include <vector>

void prep_bed(
    std::string bed_path,
    std::string bim_path,
    std::string fam_path,
    std::string out_dir,
    size_t mem_gb);
