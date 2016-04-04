#include <iostream>
#include <fstream>
#include <string>
#include <string.h>
using namespace std;
 
void output(string& line) {
	char buf[line.size() + 1];
	strcpy(buf, line.c_str());
	char* start = buf;
	char* l = buf;
	int count = 0;
	while (*l) {
		if (*l == ',') {
			count++;
			*l = 0;
			if (count == 2 || count == 3 || count == 8)
				cout << start;
			else
				cout << '"' << start << '"';
			cout << ',';
			start = l + 1;
		}
		l++;
	}
	cout << '"' << start << "\"\n";
}

int main(int argc, char* argv[]) {
	char* chr = argv[2];
	int chr_len = strlen(chr);
	if (!strcmp(chr, "Any"))
		chr = NULL;

	char snp_type = argv[3][0];
	char from = 0;
	char to;
	if (snp_type != 'd' && snp_type != 'i') {
		from = snp_type;
		to = argv[3][2];
		snp_type = 's';
	}

	string line;
	ifstream f (argv[1]);
	while (getline(f, line)) {
		if (chr)
			if (strncmp(line.c_str(), chr, chr_len) || line[chr_len] != ',')
				continue;
		int i = 0;
		for (int n=line.size(),c=0; i<n && c<3; i++) {
			if (line[i] == ',')
				++c;
		}
		if (line[i] == snp_type)
			if (from) {
				int comma = line.find(',', i + 4);
				if (from == line[comma + 1] && to == line[comma + 3])
					output(line);
			} else
				output(line);
	}
	f.close();
	return 0;
}