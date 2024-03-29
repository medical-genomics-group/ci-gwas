#pragma once

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

const int BED_PREFIX_BYTES = 3;
const unsigned char BED_PREFIX_COL_MAJ[3] = {0x6c, 0x1b, 0x01};

class BfilesBase
{
   private:
    std::string base;

   public:
    BfilesBase(std::string path) : base(path) {}

    std::string dim() const { return base + ".dim"; }

    std::string bed() const { return base + ".bed"; }

    std::string means() const { return base + ".means"; }

    std::string stds() const { return base + ".stds"; }

    std::string bim() const { return base + ".bim"; }

    std::string fam() const { return base + ".fam"; }

    std::string modes() const { return base + ".modes"; }

    std::string blocks() const { return base + ".blocks"; }

    std::string blocks(int size) const { return base + "_m" + std::to_string(size) + ".blocks"; }

    std::string with_suffix(std::string sfx) const { return base + "." + sfx; }

    bool has_valid_bed_prefix() const
    {
        std::vector<unsigned char> buffer(BED_PREFIX_BYTES);
        std::ifstream bed_file(bed(), std::ios::binary);
        bed_file.read(reinterpret_cast<char *>(buffer.data()), BED_PREFIX_BYTES);
        for (size_t i = 0; i < BED_PREFIX_BYTES; ++i)
        {
            if (buffer[i] != BED_PREFIX_COL_MAJ[i])
            {
                return false;
            }
        }
        return true;
    }
};
