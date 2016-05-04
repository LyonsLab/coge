#include <fstream>
#include <iostream>
#include <string>
#include <string.h>
#include <vector>
using namespace std;

int nth_token(string& s, int n) {
	int i = 0;
	int total_commas = n - 1;
	for (int num_commas=0,l=s.size(); i<l && num_commas<total_commas; i++) {
		if (s[i] == ',')
			num_commas++;
	}
	return i;
}

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

void split(char* s, vector<string>& tokens) {
	char* p = s + 1;
	while (char c = *p) {
		if (c == ',') {
			tokens.push_back(string(s, p - s));
			s = p + 1;
			p = s + 1;
		} else
			p++;
	}
	tokens.push_back(string(s));
}

// usage: snp_search file chromosome snp_types
// snp_types is comma delimited
int main(int argc, char* argv[]) {
	char* chr = argv[2];
	int chr_len = strlen(chr);
	if (!strcmp(chr, "Any"))
		chr = NULL;

	vector<string> snp_types;
	split(argv[3], snp_types);

	string line;
	ifstream f (argv[1]);
	while (getline(f, line)) {
		if (chr)
			if (strncmp(line.c_str(), chr, chr_len) || line[chr_len] != ',')
				continue;

		int fourth_token = nth_token(line, 4);
		for (vector<string>::iterator snp_type = snp_types.begin(); snp_type != snp_types.end(); snp_type++) {
			char initial = (*snp_type)[0];
			char from = 0;
			char to;
			if (initial != 'd' && initial != 'i') {
				from = initial;
				to = (*snp_type)[2];
				initial = 's';
			}
			if (line[fourth_token] == initial)
				if (from) {
					int comma = line.find(',', fourth_token + 4);
					if (from == line[comma + 1] && to == line[comma + 3])
						output(line);
				} else
					output(line);
		}
	}
	f.close();
	return 0;
}
