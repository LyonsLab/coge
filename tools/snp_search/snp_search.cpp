#include <iostream>
#include <fstream>
#include <string>
using namespace std;
 
 
int main(int argc, char* argv[]) {
	char snp_type = argv[2][0];
	char from = 0;
	char to;
	if (snp_type != 'd' && snp_type != 'i') {
		from = snp_type;
		to = argv[2][2];
		snp_type = 's';
	}

	string line;
	ifstream f (argv[1]);
	while (getline(f, line)) {
		int i = 0;
		for (int n=line.size(),c=0; i<n && c<3; i++) {
			if (line[i] == ',')
				++c;
		}
		if (line[i] == snp_type)
			if (from) {
				if (from == line[i + 6] && to == line[i + 8])
					cout << line << '\n';
			} else
				cout << line << '\n';
	}
	f.close();
	return 0;
}