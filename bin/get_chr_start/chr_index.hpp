
#include <map>
#include <string>
using namespace std;

class chr_index {
public:
	chr_index(char* filename);

	unsigned int	get_offset(char* chromosome);
	unsigned int	get_num_lines(char* chromosome);

private:
	void	read_index();
	void	write_index();

	string						m_filename;
	string						m_index_filename;
	map<string,unsigned int>	m_num_lines;
	map<string,unsigned int>	m_offsets;
};