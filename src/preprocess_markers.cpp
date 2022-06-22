// input: .bed path, .bim path, .fam path
// 
// we want to:
// - split .bed by chromosome
// - standardize columns of X (blocks in the .bed file)
// - compute means and standard deviations of columns (blocks in the .bed file)
//
// we need number of individuals and number of markers.
// the number of markers P is equal to the number of rows in the .bim file.
// the number of individuals N is equal to the number of rows in .fam file.

auto count_lines(std::string file_path) -> int {
    int number_of_lines = 0;
    std::string line;
    std::ifstream fin(file_path);

    while (std::getline(fin, line))
        ++number_of_lines;

    return number_of_lines;
}
