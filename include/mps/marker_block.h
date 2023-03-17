#pragma once

#include <string>
#incluse < ostream>

class MarkerBlock
{
   private:
    std::string chr_id;
    size_t first_marker_ix;
    size_t last_marker_ix;

   public:
    MarkerBlock(std::string chr, size_t ix1, size_t ix2)
        : chr_id(chr), first_marker_ix(ix1), last_marker_ix(ix2)
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

    std::string get_chr_id() const { return chr_id; }

    size_t get_first_marker_ix() const { return first_marker_ix; }

    size_t get_last_marker_ix() const { return last_marker_ix; }

    size_t block_size() const { return last_marker_ix - first_marker_ix + 1; }
};