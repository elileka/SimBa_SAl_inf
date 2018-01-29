#ifndef ___COMPUTE_ALIGNMENT_DP_H
#define ___COMPUTE_ALIGNMENT_DP_H	

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
#include "chop_prob.h"
#include "quick_jtt.h"
// code by Eli Levy Karin

#include "read_input_seqs.h"
using namespace std;

class compute_alignment_dp
{
public:

	compute_alignment_dp(read_chop_tables & chop_tables_obj, read_input_seqs & read_input_seqs_obj, const double branch_length_t, bool is_jc, quick_jtt & quick_jtt_obj) :
		_chop_tables_obj(chop_tables_obj), _read_input_seqs_obj(read_input_seqs_obj), _branch_length_t(branch_length_t), _is_jc(is_jc), _quick_jtt_obj(quick_jtt_obj)
	{
		_coded_seqs.push_back(read_input_seqs_obj.get_coded_anc_unaligned());
		_coded_seqs.push_back(read_input_seqs_obj.get_coded_des_unaligned());
		_length_of_anc = _coded_seqs[0].size();
		_length_of_des = _coded_seqs[1].size();	
	}

	vector<string> get_sampled_alignment(double & sampled_alignment_log_prob, bool should_take_max, int band_width); // If band_width == -1, no corner cutting (all pairs). Else, cuts using a fixed-width band
	double compute_conditional_alignment_total_log_prob(int band_width); // DP + If band_width == -1, no corner cutting (all pairs). calls fill_P_table_corner_cutting() and then - eq 15 of Miklos et al 2004 with zero indices
	vector<vector<double>> get_P_table(int band_width) { if (!_is_P_computed) { fill_P_table_corner_cutting(band_width); } return _P_table; } // DP + If band_width == -1, no corner cutting (all pairs). Else, cuts using a fixed-width band
	vector<vector<double>> get_X_table(int band_width) { if (!_is_X_computed) { fill_X_table_corner_cutting(band_width); } return _X_table; } // DP + If band_width == -1, no corner cutting (all pairs). Else, cuts using a fixed-width band


private:
	const double _branch_length_t;
	read_chop_tables & _chop_tables_obj;
	quick_jtt & _quick_jtt_obj;
	read_input_seqs & _read_input_seqs_obj;

	vector<vector<size_t>> _coded_seqs; // 1 - 'A', 2 - 'C', 3 - 'G', 4 - 'T' for DNA, similarly for AA
	size_t _length_of_anc; // number of chars in ancestor (first) seq (unaligned)
	size_t _length_of_des; // number of chars in descendant (second) seq  (unaligned)
	bool _is_jc;

	bool _is_P_computed = false; // to avoid recomputing if called a second time
	bool _is_X_computed = false; // to avoid recomputing if called a second time

	double _conditional_alignment_log_total_prob = 0.0; // log(Pr(B|A)) based on eq 14 + eq 15 of Miklos et al 2004

	vector<vector<double>> _P_table; // contains the "forward" probabilities (eq 14 of Miklos et al 2004 with zero indices)
	void fill_P_table_corner_cutting(int band_width); // DP + If band_width == -1, no corner cutting (all pairs). Else, cuts using a fixed-width band
	
	vector<vector<double>> _X_table; // contains the "backward" probabilities (analog of eq 14 of Miklos et al 2004 with zero indices)
	void fill_X_table_corner_cutting(int band_width); // DP + If band_width == -1, no corner cutting (all pairs). Else, cuts using a fixed-width band

	vector<vector<pair<chop_prob, double>>> _S_table; // a sampled path based on a modification of eq 14 of Miklos et al 2004 with zero indices
	void fill_S_table_corner_cutting(bool should_take_max, int band_width); // DP + corener cutting (if band_width == -1, no corner cutting); if should_take_max == true, the sample will be of the most probable
	chop_prob get_sampled_path_last_chop(bool should_take_max, double & sampled_alignment_log_prob); // based on eq 15 of Miklos et al 2004 with zero indices; if should_take_max == true, the sample will be of the most probable

	size_t sample_index_gumbel_max(vector<double> & log_probs);

	void get_combs(int band_width, vector<vector<size_t>> & pairs_to_consider); // creates a data structure of pairs to go over. If band_width == -1, no corner cutting (all pairs). Else, cuts using a fixed-width band
	void get_prev_combs_for_pair(int band_width, vector<vector<size_t>> & prev_pairs_to_consider, size_t anc_ind_i, size_t des_ind_j); // creates a data structure of prev pairs to go over. If band_width == -1, no corner cutting (all pairs). Else, cuts using a fixed-width band
	void get_backward_combs(int band_width, vector<vector<size_t>> & pairs_to_consider_backward); // creates a data structure of pairs to go over. If band_width == -1, no corner cutting (all pairs). Else, cuts using a fixed-width band
	void get_backward_prev_combs_for_pair(int band_width, vector<vector<size_t>> & prev_pairs_to_consider_backward, size_t anc_ind_i, size_t des_ind_j); // creates a data structure of prev pairs to go over. If band_width == -1, no corner cutting (all pairs). Else, cuts using a fixed-width band

};

#endif


