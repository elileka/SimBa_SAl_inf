// code by Eli Levy Karin

#ifndef ___ALIGNMENT_H
#define ___ALIGNMENT_H	

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
#include "chop_prob.h"
#include "read_input_seqs.h"
using namespace std;

class alignment
{
public:

	alignment(read_chop_tables & chop_tables_obj, vector<vector<size_t>> & coded_alignment, const double branch_length_t, bool is_jc, quick_jtt & quick_jtt_obj) :
		_chop_tables_obj(chop_tables_obj), _coded_alignment(coded_alignment), _branch_length_t(branch_length_t), _is_jc(is_jc), _quick_jtt_obj(quick_jtt_obj)
	{
		_length_of_aligment = _coded_alignment.size();
		break_alignment_to_chops();
	} //constructor

	double get_alignment_log_probability_cond_on_anc() const;


private:
	const double _branch_length_t;
	read_chop_tables & _chop_tables_obj;
	quick_jtt & _quick_jtt_obj;
	bool _is_jc;
	//read_input_seqs & _read_input_seqs_obj;
	
	vector<vector<size_t>> _coded_alignment; // '0' denotes a gap, 1 - 'A', 2 - 'C', 3 - 'C', 4 - 'T' (for DNA)
	size_t _length_of_aligment; // number of chars in each seq

	vector<chop_prob> _aligment_chops;

	void break_alignment_to_chops();

};

#endif


