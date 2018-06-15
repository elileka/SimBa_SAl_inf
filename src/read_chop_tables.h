// code by Eli Levy Karin

#ifndef ___READ_CHOP_TABLES_H
#define ___READ_CHOP_TABLES_H	

#include <cstdio>
#include <cmath>
#include <vector>
#include <random>
#include <iostream>
#include <sstream>
#include <map>
#include <fstream>
#include <string>
using namespace std;

class read_chop_tables
{
public:

	read_chop_tables(const string chop_tables_file_string) :
		_chop_tables_file_string(chop_tables_file_string)
	{
		read_tables_from_file();
	} //constructor

	double get_chop_prob(const char chop_type, size_t i, size_t j);
	size_t get_max_i_in_L_table() const { return _max_i_in_L_table; }
	size_t get_max_j_in_L_table() const { return _max_j_in_L_table; }
	size_t get_max_i_in_R_table() const { return _max_i_in_R_table; }
	size_t get_max_j_in_R_table() const { return _max_j_in_R_table; }
	map<pair<size_t, size_t>, double> get_N_chops_map() const { return _N_chop_probs_map; }
	


private:
	const string _chop_tables_file_string;
	double _min_probability = 1.0;

	map<pair<size_t, size_t>, double> _N_chop_probs_map;
	map<pair<size_t, size_t>, double> _L_chop_probs_map;
	map<pair<size_t, size_t>, double> _R_chop_probs_map;
	map<pair<size_t, size_t>, double> _B_chop_probs_map;

	size_t _max_i_in_L_table = 0;
	size_t _max_j_in_L_table = 0;
	size_t _max_i_in_R_table = 0;
	size_t _max_j_in_R_table = 0;

	void read_tables_from_file();

};

#endif

