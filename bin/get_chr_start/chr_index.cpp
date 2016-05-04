
#include "chr_index.hpp"

#include <string.h>
#include <sys/stat.h>

#include <fstream>

//----------------------------------------

inline bool file_exists(string& file) {
    struct stat buf;
    return stat(file.c_str(), &buf) == 0;
}

//----------------------------------------

chr_index::chr_index(char* filename) : m_filename(filename), m_index_filename(filename) {
	m_index_filename += ".chr";
	if (file_exists(m_index_filename))
		read_index();
	else
		write_index();
}

//----------------------------------------

unsigned int
chr_index::get_offset(char* chromosome) {
	return m_offsets.find(string(chromosome))->second;
}

//----------------------------------------

unsigned int
chr_index::get_num_lines(char* chromosome)
{
	return m_num_lines.find(string(chromosome))->second;
}

//----------------------------------------

void
chr_index::read_index() {
	ifstream index_file(m_index_filename.c_str());
	string chr;
	unsigned int offset;
	unsigned int num_lines;
	while (index_file >> chr >> offset >> num_lines) {
		m_offsets[chr] = offset;
		m_num_lines[chr] = num_lines;
	}
	index_file.close();
}

//----------------------------------------

void
chr_index::write_index() {
	string line;
	ifstream txt_file(m_filename.c_str());
	getline(txt_file, line);
	string chr = line.substr(0, line.find(','));
	unsigned int num_lines = 0;
	unsigned int offset = 0;
	m_offsets[chr] = offset;
	while (getline(txt_file, line)) {
		if (strncmp(chr.c_str(), line.c_str(), chr.size()) || line[chr.size()] != ',') {
			m_num_lines[chr] = num_lines;
			num_lines = 0;
			chr = line.substr(0, line.find(','));
			m_offsets[chr] = offset;
		}
		offset += line.size() + 1;
		num_lines++;
	}
	m_num_lines[chr] = num_lines;
	txt_file.close();

	ofstream index_file(m_index_filename.c_str());
	for (map<string,unsigned int>::iterator i=m_offsets.begin(); i!=m_offsets.end(); i++)
		index_file << i->first << ' ' << i->second << ' ' << m_num_lines.find(i->first)->second << '\n';
	index_file.close();
}