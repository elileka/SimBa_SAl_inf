// code by Eli Levy Karin

#ifndef ___CHOP_PROB_H
#define ___CHOP_PROB_H	

#include <cstdio>
#include <cmath>
#include <vector>
#include <random>
#include <iostream>
#include <sstream>
#include <map>
#include <fstream>
#include <string>
#include "read_chop_tables.h"
#include "quick_jtt.h"
using namespace std;

class chop_prob
{
public:

	chop_prob(read_chop_tables & chop_tables_obj, char chop_type, size_t i, size_t j, size_t anc_matched_char, size_t des_matched_char, vector<size_t> & des_inserted_chars, double branch_length_t, bool is_jc, quick_jtt & quick_jtt_obj) :
		_chop_tables_obj(chop_tables_obj), _branch_length_t(branch_length_t), _chop_type(chop_type), _i(i), _j(j), _anc_matched_char(anc_matched_char), _des_matched_char(des_matched_char), _des_inserted_chars(des_inserted_chars), _is_jc(is_jc), _quick_jtt_obj(quick_jtt_obj)
	{
		compute_chop_prob();
	} //constructor
	chop_prob& operator = (const chop_prob&);

	double get_chop_log_prob() const { return _chop_log_prob; }
	void print_chop_info() const;
	char get_chop_type() const { return _chop_type; }
	size_t get_i() const { return _i; }
	size_t get_j() const { return _j; }

private:
	read_chop_tables & _chop_tables_obj;
	quick_jtt & _quick_jtt_obj;
	double _branch_length_t;
	char _chop_type;
	bool _is_jc;

	size_t _i;
	size_t _j;
	
	vector<size_t> _des_inserted_chars; // this will be a vector of j chars

	size_t _anc_matched_char; // if chop is of type L or N, it will be different than 0
	size_t _des_matched_char; // if chop is of type L or N, it will be different than 0
	
	double _chop_log_prob = 0; // _sub_log_prob + _indel_log_prob

	double _sub_log_prob = 0;
	void compute_JC_substitution_prob();

	void compute_JTT_substitution_prob();

	double _indel_log_prob = 0;
	void compute_chop_prob();

};

#endif

