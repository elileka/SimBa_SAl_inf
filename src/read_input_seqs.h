// code by Eli Levy Karin

#ifndef ___READ_INPUT_SEQS_H
#define ___READ_INPUT_SEQS_H	

#include <cstdio>
#include <cmath>
#include <vector>
#include <random>
#include <iostream>
#include <sstream>
#include <map>
#include <fstream>
#include <string>
#include <queue>
#include <stack>
#include <climits>
using namespace std;

class read_input_seqs
{
public:

	read_input_seqs(const string seqs_fasta_file, bool is_aligned, bool is_dna) :
		_is_aligned(is_aligned), _is_dna(is_dna)
	{
		read_seqs_from_fasta_file(seqs_fasta_file);
		get_original_unaligned_seqs();
		code_seqs(_original_seqs, _coded_seqs);
		get_coded_unaligned_seqs();
		if (is_aligned)
		{
			transpose_to_get_alignment(_coded_seqs, _coded_alignment);
		}
	} //constructor

	string get_str_orig_seq(size_t seq_ind) { return _original_seqs[seq_ind]; };
	vector<vector<size_t>> get_coded_alignment() { return _coded_alignment; };
	size_t get_length_of_alignment() { return _length_of_anc; };

	string get_str_from_coded_alignment(size_t seq_ind, vector<vector<size_t>> coded_alignment);
	size_t get_length_of_anc_unaligned() { return _length_of_anc_unaligned; };
	size_t get_length_of_des_unaligned() { return _length_of_des_unaligned; };

	vector<size_t> get_coded_anc_unaligned() { return _coded_seqs_unaligned[0]; };
	vector<size_t> get_coded_des_unaligned() { return _coded_seqs_unaligned[1]; };

	void compute_gotoh_alignment();
	vector<vector<size_t>> get_gotoh_coded_alignment() { return _gotoh_coded_alignment; };
	double get_estimated_branch_length_MP() { return _branch_length_estimate_MP; };

private:
	bool _is_aligned;
	bool _is_dna;

	map<string, string> _original_seqs_and_names; // name to seq, DNA chars + "-" or AA chars. Seqs can be unaligned or aligned
	vector<string> _original_seqs; // sequences only, no headers. Seqs can be unaligned or aligned
	vector<vector<size_t>> _coded_seqs; // '0' denotes a gap, 1 - 'A', 2 - 'C', 3 - 'C', 4 - 'T' (for DNA)
	size_t _length_of_anc; // length of anc unaligned / aligned
	size_t _length_of_des; // length of des unaligned / aligned
	
	vector<string> _original_seqs_unaligned; // sequences only, no headers. Seqs are unaligned
	vector<vector<size_t>> _coded_seqs_unaligned; // 1 - 'A', 2 - 'C', 3 - 'C', 4 - 'T' (for DNA)
	size_t _length_of_anc_unaligned; // length of anc unaligned
	size_t _length_of_des_unaligned; // length of des unaligned

	// if the input is aligned, we also transpose its code: instead of item 0 being the ancestral sequence vector
	// item 0 will be position 0 in the alignment (a vector of anc and des)
	vector<vector<size_t>> _coded_alignment; // '0' denotes a gap, 1 - 'A', 2 - 'C', 3 - 'C', 4 - 'T' (for DNA)

	void read_seqs_from_fasta_file(const string seqs_fasta_file);
	void get_original_unaligned_seqs();
	void code_seqs(vector<string> & seqs_to_code, vector<vector<size_t>> & seqs_after_code);
	void get_coded_unaligned_seqs();
	void transpose_to_get_alignment(vector<vector<size_t>> & coded_seqs_before_transpose, vector<vector<size_t>> & coded_alignment_after_transpose);
	
	vector<string> _gotoh_aligned_seqs;
	vector<vector<size_t>> _gotoh_coded_seqs;
	vector<vector<size_t>> _gotoh_coded_alignment; // '0' denotes a gap, 1 - 'A', 2 - 'C', 3 - 'C', 4 - 'T' (for DNA)
	
	int _default_match_score = 2;
	int _default_mismatch_score = -1;
	int _gap_open = 5;
	int _gap_extend = 1;
	double _branch_length_estimate_MP = 0;
	map<pair<char, char>, int> _similarityMap;
	vector<string> computeAffineNWForPair(const string & A, const string & B);
	void initNucMap();
	void initAAMap();
	int get_pair_score(char a, char b);
	void estimate_branch_MP();
	
};

#endif

