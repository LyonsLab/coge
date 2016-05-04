
#include "chr_index.hpp"

#include <iostream>

// usage: get_chr_start file chromosome
int main(int argc, char* argv[]) {
	chr_index i(argv[1]);
	cout << i.get_offset(argv[2]) << ' ' << i.get_num_lines(argv[2]);
	return 0;
}