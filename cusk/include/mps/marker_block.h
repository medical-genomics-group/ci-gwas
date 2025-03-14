#pragma once

#include <ostream>
#include <sstream>
#include <string>

class MarkerBlock
{
   private:
    std::string chr_id;
    // index on the chromosome
    size_t first_marker_ix;
    size_t last_marker_ix;
    size_t chr_global_offset;

   public:
    MarkerBlock() {}

    MarkerBlock(std::string chr, size_t ix1, size_t ix2, size_t offset)
        : chr_id(chr), first_marker_ix(ix1), last_marker_ix(ix2), chr_global_offset(offset)
    {
    }

    friend std::ostream& operator<<(std::ostream& os, const MarkerBlock& bar)
    {
        return os << "MarkerBlock(" << bar.chr_id << ", " << bar.first_marker_ix << ", "
                  << bar.last_marker_ix << ")";
    }

    bool operator==(const MarkerBlock& other) const
    {
        return (this->chr_id == other.chr_id) && (this->first_marker_ix == other.first_marker_ix) &&
               (this->last_marker_ix == other.last_marker_ix);
    }

    std::string to_line_string() const
    {
        std::ostringstream ss;
        ss << chr_id << "\t" << first_marker_ix << "\t" << last_marker_ix;
        return ss.str();
    }

    std::string to_file_string() const
    {
        std::ostringstream ss;
        ss << chr_id << "_" << first_marker_ix << "_" << last_marker_ix;
        return ss.str();
    }

    std::string get_chr_id() const { return chr_id; }

    size_t get_first_marker_ix() const { return first_marker_ix; }

    size_t get_last_marker_ix() const { return last_marker_ix; }

    size_t get_first_marker_global_ix() const { return first_marker_ix + chr_global_offset; }

    size_t get_last_marker_global_ix() const { return last_marker_ix + chr_global_offset; }

    size_t block_size() const { return last_marker_ix - first_marker_ix + 1; }
};
