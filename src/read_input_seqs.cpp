// code by Eli Levy Karin

#include "read_input_seqs.h"

string read_input_seqs::get_str_from_coded_alignment(size_t seq_ind, vector<vector<size_t>> coded_alignment)
{
	string aligned_seq_str;
	for (size_t alignment_ind = 0; alignment_ind < coded_alignment.size(); alignment_ind++)
	{
		size_t seq_char = coded_alignment[alignment_ind][seq_ind];
		if (_is_dna)
		{
			switch (seq_char)
			{
			case 0: aligned_seq_str.append("-"); break;
			case 1: aligned_seq_str.append("A"); break;
			case 2: aligned_seq_str.append("C"); break;
			case 3: aligned_seq_str.append("G"); break;
			case 4: aligned_seq_str.append("T"); break;
			}
		}
		else
		{
			switch (seq_char)
			{
			case 0: aligned_seq_str.append("-"); break;
			case 1: aligned_seq_str.append("A"); break;
			case 2: aligned_seq_str.append("R"); break;
			case 3: aligned_seq_str.append("N"); break;
			case 4: aligned_seq_str.append("D"); break;
			case 5: aligned_seq_str.append("C"); break;
			case 6: aligned_seq_str.append("Q"); break;
			case 7: aligned_seq_str.append("E"); break;
			case 8: aligned_seq_str.append("G"); break;
			case 9: aligned_seq_str.append("H"); break;
			case 10: aligned_seq_str.append("I"); break;
			case 11: aligned_seq_str.append("L"); break;
			case 12: aligned_seq_str.append("K"); break;
			case 13: aligned_seq_str.append("M"); break;
			case 14: aligned_seq_str.append("F"); break;
			case 15: aligned_seq_str.append("P"); break;
			case 16: aligned_seq_str.append("S"); break;
			case 17: aligned_seq_str.append("T"); break;
			case 18: aligned_seq_str.append("W"); break;
			case 19: aligned_seq_str.append("Y"); break;
			case 20: aligned_seq_str.append("V"); break;
			}
		}

	}
	return aligned_seq_str;
}

void read_input_seqs::read_seqs_from_fasta_file(const string seqs_fasta_file)
{
	ifstream myFile;
	myFile.open(seqs_fasta_file.c_str());
	if (!myFile)
	{
		cerr << "can not open file: " << seqs_fasta_file << endl;
		exit(1);
	}

	string name = "";
	string seq = "";

	string fasta_line;
	while (getline(myFile, fasta_line))
	{
		//remove white spaces from end of 'fasta_line':
		fasta_line.erase(fasta_line.find_last_not_of(" \n\r\t") + 1);

		if ((fasta_line[0] == '>') || (name.size() > 1))
		{
			if (fasta_line[0] == '>') // if this is the first sequence
			{
				name = fasta_line.substr(1);
			}
			else // this is not the first sequence
			{
				seq += fasta_line;
			}
			while (getline(myFile, fasta_line))
			{
				//remove white spaces from end of 'fasta_line':
				fasta_line.erase(fasta_line.find_last_not_of(" \n\r\t") + 1);

				if (fasta_line[0] == '>') // reached the next seq
				{
					_original_seqs_and_names[name] = seq; // add what you have
					_original_seqs.push_back(seq); // collect just the aligned sequence - not the header
					name = fasta_line.substr(1); // update name
					seq = ""; // initialize seq
					break;
				}
				else
				{
					seq += fasta_line;
				}
			}
		}
	}
	_original_seqs_and_names[name] = seq; // insert the last sequence
	_original_seqs.push_back(seq); // collect just the aligned sequence - not the header
	myFile.close();

	// validate there are two sequences:
	size_t num_seqs = _original_seqs.size();
	if (num_seqs != 2)
	{
		cerr << "more than two sequences were provided: " << seqs_fasta_file << endl;
		exit(1);
	}
	
	// validate all sequences are longer than 0. if aligned - same length!
	size_t length_of_first_seq = _original_seqs[0].size();
	for (size_t i = 0; i < _original_seqs.size(); i++)
	{
		if (_original_seqs[i].size() == 0)
		{
			cerr << "at least one sequence is of length 0 - we don't handle this: " << seqs_fasta_file << endl;
			exit(1);
		}
		if (_is_aligned)
		{
			if (_original_seqs[i].size() != length_of_first_seq)
			{
				cerr << "sequences are not all of the same length: " << seqs_fasta_file << endl;
				exit(1);
			}
		}
	}
	
	_length_of_anc = _original_seqs[0].size();
	_length_of_des = _original_seqs[1].size();
}

void read_input_seqs::get_original_unaligned_seqs()
{
	if (!_is_aligned)
	{
		_original_seqs_unaligned.push_back(_original_seqs[0]);
		_original_seqs_unaligned.push_back(_original_seqs[1]);
	}
	else
	{
		string anc_seq_unaligned = "";
		string des_seq_unaligned = "";
		for (char& c : _original_seqs[0])
		{
			if (c != '-')
			{
				anc_seq_unaligned = anc_seq_unaligned + c;
			}
		}
		for (char& c : _original_seqs[1])
		{
			if (c != '-')
			{
				des_seq_unaligned = des_seq_unaligned + c;
			}
		}
		_original_seqs_unaligned.push_back(anc_seq_unaligned);
		_original_seqs_unaligned.push_back(des_seq_unaligned);
	}
}

void read_input_seqs::code_seqs(vector<string> & seqs_to_code, vector<vector<size_t>> & seqs_after_code)
{
	for (size_t i = 0; i < seqs_to_code.size(); i++)
	{
		vector<size_t> curr_coded_seq;
		string sequence = seqs_to_code[i];

		for (char & seq_char : sequence)
		{
			size_t char_code;
			if (_is_dna)
			{
				switch (seq_char)
				{
				case '-': char_code = 0; break;
				case 'A': case'a': char_code = 1; break;
				case 'C': case'c': char_code = 2; break;
				case 'G': case'g': char_code = 3; break;
				case 'T': case't': char_code = 4; break;
				default:
					cerr << "illegal char: " << seq_char << ", should be [AaCcGgTt-]" << endl;
					exit(1);
				}
			}
			else
			{
				// amino acids
				switch (seq_char)
				{
				case '-': char_code = 0; break;
				case 'A': case'a': char_code = 1; break;
				case 'R': case'r': char_code = 2; break;
				case 'N': case'n': char_code = 3; break;
				case 'D': case'd': char_code = 4; break;
				case 'C': case'c': char_code = 5; break;
				case 'Q': case'q': char_code = 6; break;
				case 'E': case'e': char_code = 7; break;
				case 'G': case'g': char_code = 8; break;
				case 'H': case'h': char_code = 9; break;
				case 'I': case'i': char_code = 10; break;
				case 'L': case'l': char_code = 11; break;
				case 'K': case'k': char_code = 12; break;
				case 'M': case'm': char_code = 13; break;
				case 'F': case'f': char_code = 14; break;
				case 'P': case'p': char_code = 15; break;
				case 'S': case's': char_code = 16; break;
				case 'T': case't': char_code = 17; break;
				case 'W': case'w': char_code = 18; break;
				case 'Y': case'y': char_code = 19; break;
				case 'V': case'v': char_code = 20; break;
				default:
					cerr << "illegal char: " << seq_char << ", should be [AaRrNnDdCcQqEeGgHhIiLlKkMmFfPpSsTtWwYyVv-]" << endl;
					exit(1);
				}
			}
			curr_coded_seq.push_back(char_code);
		}
		seqs_after_code.push_back(curr_coded_seq);

	}
	
}

void read_input_seqs::get_coded_unaligned_seqs()
{
	if (!_is_aligned)
	{
		_coded_seqs_unaligned.push_back(_coded_seqs[0]);
		_coded_seqs_unaligned.push_back(_coded_seqs[1]);
	}
	else
	{
		vector<size_t> anc_coded_seq_unaligned;
		vector<size_t> des_coded_seq_unaligned;
		size_t length_of_alignment = _coded_seqs[0].size();

		for (size_t alignment_ind = 0; alignment_ind < length_of_alignment; alignment_ind++)
		{
			size_t anc_code = _coded_seqs[0][alignment_ind];
			size_t des_code = _coded_seqs[1][alignment_ind];

			if (anc_code != 0) // 0 is a gap
			{
				anc_coded_seq_unaligned.push_back(anc_code);
			}
			if (des_code != 0) // 0 is a gap
			{
				des_coded_seq_unaligned.push_back(des_code);
			}
		}
		_coded_seqs_unaligned.push_back(anc_coded_seq_unaligned);
		_coded_seqs_unaligned.push_back(des_coded_seq_unaligned);

		_length_of_anc_unaligned = anc_coded_seq_unaligned.size();
		_length_of_des_unaligned = des_coded_seq_unaligned.size();

	}
}

void read_input_seqs::transpose_to_get_alignment(vector<vector<size_t>> & coded_seqs_before_transpose, vector<vector<size_t>> & coded_alignment_after_transpose)
{
	// now let us transpose the alignment:
	size_t length_of_aligment = coded_seqs_before_transpose[0].size();
	for (size_t pos_ind = 0; pos_ind < length_of_aligment; pos_ind++)
	{
		vector<size_t> curr_pos;
		for (size_t seq_ind = 0; seq_ind < coded_seqs_before_transpose.size(); seq_ind++)
		{
			curr_pos.push_back(coded_seqs_before_transpose[seq_ind][pos_ind]);
		}
		coded_alignment_after_transpose.push_back(curr_pos);
	}
}

void read_input_seqs::compute_gotoh_alignment()
{
	string anc_unaligned = _original_seqs_unaligned[0];
	string des_unaligned = _original_seqs_unaligned[1];

	if (_is_dna)
	{
		initNucMap();
	}
	else
	{
		initAAMap();
	}

	_gotoh_aligned_seqs = computeAffineNWForPair(anc_unaligned, des_unaligned);
	code_seqs(_gotoh_aligned_seqs, _gotoh_coded_seqs);
	transpose_to_get_alignment(_gotoh_coded_seqs, _gotoh_coded_alignment);
	estimate_branch_MP();
}

void read_input_seqs::estimate_branch_MP()
{
	size_t num_non_gap_positions = 0;
	size_t num_diff_positions = 0;
	for (size_t i = 0; i < _gotoh_coded_alignment.size(); i++)
	{
		if ((_gotoh_coded_alignment[i][0] != 0) && (_gotoh_coded_alignment[i][1] != 0))
		{
			num_non_gap_positions++;
			if (_gotoh_coded_alignment[i][0] != _gotoh_coded_alignment[i][1])
			{
				num_diff_positions++;
			}
		}
	}
	if (num_non_gap_positions != 0)
	{
		_branch_length_estimate_MP = (double)num_diff_positions / (double)num_non_gap_positions;
	}
}

vector<string> read_input_seqs::computeAffineNWForPair(const string & A, const string & B)
{
	vector<string> aligned_seqs;

	size_t m = A.length();
	size_t n = B.length();

	// deal with edge cases where one or more of the sequences is of size zero:
	if ((A.length() == 0) && (B.length() == 0))
	{

		aligned_seqs.push_back("");
		aligned_seqs.push_back("");
		return aligned_seqs;
	}

	if (A.length() == 0)
	{
		string retA = "";
		retA.insert(0, B.length(), '-'); // if A is empty we enter B.length() instances of '-' to the aligned sequence.
		aligned_seqs.push_back(retA);
		aligned_seqs.push_back(B);
		//int aln_score = _gap_open + (B.length() - 1)*_gap_extend;
		return aligned_seqs;
	}

	if (B.length() == 0)
	{
		string retB = "";
		retB.insert(0, A.length(), '-');
		aligned_seqs.push_back(A);
		aligned_seqs.push_back(retB);
		//int aln_score = _gap_open + (A.length() - 1)*_gap_extend;
		return aligned_seqs;
	}
	// end deal with edge cases where one or more of the sequences is of size zero
	// This computation is based on the book by Durbin (Biological sequence analysis), pg 29 

	// The M matrix. M[i][j] is the best score when the first sequence is aligned till position i,
	// the second sequence is aligned till position j, and position i is aligned to position j.
	vector<vector<int>> M_mat(m + 1, vector<int>(n + 1));

	// The I matrix. I[i][j] is the best score when the first sequence is aligned till position i,
	// the second sequence is aligned till position j, and position i in the first sequence is
	// aligned to a gap.
	vector<vector<int>> I_mat(m + 1, vector<int>(n + 1));

	vector<vector<int>> H_mat(m + 1, vector<int>(n + 1));

	vector<vector<int>> J_mat(m + 1, vector<int>(n + 1));

	//matrices initialization

	// Initialization of the M matrix. In M, we assume that position i is aligned to position j.
	// It cannot be that position 2 of the first sequence is aligned with a match to position 0
	// of the second sequence. So, we set M[0,i] and M[j,0] to be minus infinity.

	// The I matrix captures the cases in which the upper sequence is extended
	// For example    AAKC
	// Aligned with   A_T_
	// I[0][0] is impossible and hence it is set to minus infinity.
	// I[0][5] is also impossible and hence it is set to minus inf as well.
	// I[5][0] suggests that there are five gaps in the second sequence, so it
	// is set to gap openning + the cost of four gaps.

	// The matrix H stores the max of the matrices M,I,J.
	for (size_t j = 1; j <= B.length(); j++)
	{
		M_mat[0][j] = -(INT_MAX / 2);
		I_mat[0][j] = -(INT_MAX / 2);
		J_mat[0][j] = -_gap_open - _gap_extend *(j - 1);
		H_mat[0][j] = -_gap_open - _gap_extend *(j - 1);
	}

	for (size_t i = 1; i <= m; i++)
	{
		M_mat[i][0] = -(INT_MAX / 2);
		I_mat[i][0] = -_gap_open - _gap_extend *(i - 1);
		J_mat[i][0] = -(INT_MAX / 2);
		H_mat[i][0] = -_gap_open - _gap_extend *(i - 1);
	}

	H_mat[0][0] = 0;
	M_mat[0][0] = 0; // this is a special case of M, which is ok (zero can be aligned to zero).
	I_mat[0][0] = -(INT_MAX / 2);
	J_mat[0][0] = -(INT_MAX / 2);
	// end initialization

	//score computation - dynamic programming
	for (size_t i = 1; i <= m; i++)
	{
		for (size_t j = 1; j <= n; j++)
		{
			int S = get_pair_score(A[i - 1], B[j - 1]);
			int M_diag = M_mat[i - 1][j - 1] + S;
			int I_diag = I_mat[i - 1][j - 1] + S;
			int J_diag = J_mat[i - 1][j - 1] + S;
			M_mat[i][j] = max(M_diag, max(I_diag, J_diag));

			// here you insert a gap in the top sequence
			int M_up = M_mat[i - 1][j];
			int I_up = I_mat[i - 1][j];
			I_mat[i][j] = max((M_up - _gap_open), (I_up - _gap_extend));

			int M_left = M_mat[i][j - 1];
			int J_left = J_mat[i][j - 1];
			J_mat[i][j] = max((M_left - _gap_open), (J_left - _gap_extend));

			H_mat[i][j] = max(M_mat[i][j], max(I_mat[i][j], J_mat[i][j]));
		}
	}
	// end score computation - dynamic programming

	// traceback
	string retA, retB;
	stack<char> SA, SB;

	int ii = m;
	int jj = n;

	while (ii != 0 || jj != 0)
	{
		if (ii == 0)
		{
			SA.push('-');
			SB.push(B[jj - 1]);
			jj--;
		}
		else if (jj == 0)
		{
			SA.push(A[ii - 1]);
			SB.push('-');
			ii--;
		}
		else
		{
			int curr_H_val = H_mat[ii][jj];
			if (curr_H_val == M_mat[ii][jj])
			{
				//go diag
				SA.push(A[ii - 1]); // A [ii-1] is the ii'th position in sequence A.
				SB.push(B[jj - 1]);
				ii--;
				jj--;
			}
			else if (curr_H_val == I_mat[ii][jj])
			{
				//take from A
				SA.push(A[ii - 1]);
				SB.push('-');
				ii--;
			}
			else
			{
				//take from B
				SA.push('-');
				SB.push(B[jj - 1]);
				jj--;
			}
		}
	}

	while (!SA.empty())
	{
		retA += SA.top();
		retB += SB.top();
		SA.pop();
		SB.pop();
	}
	aligned_seqs.push_back(retA);
	aligned_seqs.push_back(retB);
	//cout << "The score of the best alignment = " << H_mat[m][n] << endl;
	return aligned_seqs;
}

void read_input_seqs::initNucMap()
{
	_similarityMap[make_pair('a', 'a')] = 10;
	_similarityMap[make_pair('a', 'c')] = -10;
	_similarityMap[make_pair('a', 'g')] = -5;
	_similarityMap[make_pair('a', 't')] = -10;

	_similarityMap[make_pair('c', 'a')] = -10;
	_similarityMap[make_pair('c', 'c')] = 10;
	_similarityMap[make_pair('c', 'g')] = -10;
	_similarityMap[make_pair('c', 't')] = -5;

	_similarityMap[make_pair('g', 'a')] = -5;
	_similarityMap[make_pair('g', 'c')] = -10;
	_similarityMap[make_pair('g', 'g')] = 10;
	_similarityMap[make_pair('g', 't')] = -10;

	_similarityMap[make_pair('t', 'a')] = -10;
	_similarityMap[make_pair('t', 'c')] = -5;
	_similarityMap[make_pair('t', 'g')] = -10;
	_similarityMap[make_pair('t', 't')] = 10;

}

void read_input_seqs::initAAMap()
{
	//#  Matrix made by matblas from blosum62.iij
	// downloaded from https://www.ncbi.nlm.nih.gov/Class/FieldGuide/BLOSUM62.txt
	// prepared for hard coding by perl script readSimMatrix.pl
	_similarityMap[make_pair('a', 'a')] = 4;
	_similarityMap[make_pair('r', 'a')] = -1;
	_similarityMap[make_pair('n', 'a')] = -2;
	_similarityMap[make_pair('d', 'a')] = -2;
	_similarityMap[make_pair('c', 'a')] = 0;
	_similarityMap[make_pair('q', 'a')] = -1;
	_similarityMap[make_pair('e', 'a')] = -1;
	_similarityMap[make_pair('g', 'a')] = 0;
	_similarityMap[make_pair('h', 'a')] = -2;
	_similarityMap[make_pair('i', 'a')] = -1;
	_similarityMap[make_pair('l', 'a')] = -1;
	_similarityMap[make_pair('k', 'a')] = -1;
	_similarityMap[make_pair('m', 'a')] = -1;
	_similarityMap[make_pair('f', 'a')] = -2;
	_similarityMap[make_pair('p', 'a')] = -1;
	_similarityMap[make_pair('s', 'a')] = 1;
	_similarityMap[make_pair('t', 'a')] = 0;
	_similarityMap[make_pair('w', 'a')] = -3;
	_similarityMap[make_pair('y', 'a')] = -2;
	_similarityMap[make_pair('v', 'a')] = 0;
	_similarityMap[make_pair('b', 'a')] = -2;
	_similarityMap[make_pair('z', 'a')] = -1;
	_similarityMap[make_pair('x', 'a')] = 0;
	_similarityMap[make_pair('*', 'a')] = -4;
	_similarityMap[make_pair('a', 'r')] = -1;
	_similarityMap[make_pair('r', 'r')] = 5;
	_similarityMap[make_pair('n', 'r')] = 0;
	_similarityMap[make_pair('d', 'r')] = -2;
	_similarityMap[make_pair('c', 'r')] = -3;
	_similarityMap[make_pair('q', 'r')] = 1;
	_similarityMap[make_pair('e', 'r')] = 0;
	_similarityMap[make_pair('g', 'r')] = -2;
	_similarityMap[make_pair('h', 'r')] = 0;
	_similarityMap[make_pair('i', 'r')] = -3;
	_similarityMap[make_pair('l', 'r')] = -2;
	_similarityMap[make_pair('k', 'r')] = 2;
	_similarityMap[make_pair('m', 'r')] = -1;
	_similarityMap[make_pair('f', 'r')] = -3;
	_similarityMap[make_pair('p', 'r')] = -2;
	_similarityMap[make_pair('s', 'r')] = -1;
	_similarityMap[make_pair('t', 'r')] = -1;
	_similarityMap[make_pair('w', 'r')] = -3;
	_similarityMap[make_pair('y', 'r')] = -2;
	_similarityMap[make_pair('v', 'r')] = -3;
	_similarityMap[make_pair('b', 'r')] = -1;
	_similarityMap[make_pair('z', 'r')] = 0;
	_similarityMap[make_pair('x', 'r')] = -1;
	_similarityMap[make_pair('*', 'r')] = -4;
	_similarityMap[make_pair('a', 'n')] = -2;
	_similarityMap[make_pair('r', 'n')] = 0;
	_similarityMap[make_pair('n', 'n')] = 6;
	_similarityMap[make_pair('d', 'n')] = 1;
	_similarityMap[make_pair('c', 'n')] = -3;
	_similarityMap[make_pair('q', 'n')] = 0;
	_similarityMap[make_pair('e', 'n')] = 0;
	_similarityMap[make_pair('g', 'n')] = 0;
	_similarityMap[make_pair('h', 'n')] = 1;
	_similarityMap[make_pair('i', 'n')] = -3;
	_similarityMap[make_pair('l', 'n')] = -3;
	_similarityMap[make_pair('k', 'n')] = 0;
	_similarityMap[make_pair('m', 'n')] = -2;
	_similarityMap[make_pair('f', 'n')] = -3;
	_similarityMap[make_pair('p', 'n')] = -2;
	_similarityMap[make_pair('s', 'n')] = 1;
	_similarityMap[make_pair('t', 'n')] = 0;
	_similarityMap[make_pair('w', 'n')] = -4;
	_similarityMap[make_pair('y', 'n')] = -2;
	_similarityMap[make_pair('v', 'n')] = -3;
	_similarityMap[make_pair('b', 'n')] = 3;
	_similarityMap[make_pair('z', 'n')] = 0;
	_similarityMap[make_pair('x', 'n')] = -1;
	_similarityMap[make_pair('*', 'n')] = -4;
	_similarityMap[make_pair('a', 'd')] = -2;
	_similarityMap[make_pair('r', 'd')] = -2;
	_similarityMap[make_pair('n', 'd')] = 1;
	_similarityMap[make_pair('d', 'd')] = 6;
	_similarityMap[make_pair('c', 'd')] = -3;
	_similarityMap[make_pair('q', 'd')] = 0;
	_similarityMap[make_pair('e', 'd')] = 2;
	_similarityMap[make_pair('g', 'd')] = -1;
	_similarityMap[make_pair('h', 'd')] = -1;
	_similarityMap[make_pair('i', 'd')] = -3;
	_similarityMap[make_pair('l', 'd')] = -4;
	_similarityMap[make_pair('k', 'd')] = -1;
	_similarityMap[make_pair('m', 'd')] = -3;
	_similarityMap[make_pair('f', 'd')] = -3;
	_similarityMap[make_pair('p', 'd')] = -1;
	_similarityMap[make_pair('s', 'd')] = 0;
	_similarityMap[make_pair('t', 'd')] = -1;
	_similarityMap[make_pair('w', 'd')] = -4;
	_similarityMap[make_pair('y', 'd')] = -3;
	_similarityMap[make_pair('v', 'd')] = -3;
	_similarityMap[make_pair('b', 'd')] = 4;
	_similarityMap[make_pair('z', 'd')] = 1;
	_similarityMap[make_pair('x', 'd')] = -1;
	_similarityMap[make_pair('*', 'd')] = -4;
	_similarityMap[make_pair('a', 'c')] = 0;
	_similarityMap[make_pair('r', 'c')] = -3;
	_similarityMap[make_pair('n', 'c')] = -3;
	_similarityMap[make_pair('d', 'c')] = -3;
	_similarityMap[make_pair('c', 'c')] = 9;
	_similarityMap[make_pair('q', 'c')] = -3;
	_similarityMap[make_pair('e', 'c')] = -4;
	_similarityMap[make_pair('g', 'c')] = -3;
	_similarityMap[make_pair('h', 'c')] = -3;
	_similarityMap[make_pair('i', 'c')] = -1;
	_similarityMap[make_pair('l', 'c')] = -1;
	_similarityMap[make_pair('k', 'c')] = -3;
	_similarityMap[make_pair('m', 'c')] = -1;
	_similarityMap[make_pair('f', 'c')] = -2;
	_similarityMap[make_pair('p', 'c')] = -3;
	_similarityMap[make_pair('s', 'c')] = -1;
	_similarityMap[make_pair('t', 'c')] = -1;
	_similarityMap[make_pair('w', 'c')] = -2;
	_similarityMap[make_pair('y', 'c')] = -2;
	_similarityMap[make_pair('v', 'c')] = -1;
	_similarityMap[make_pair('b', 'c')] = -3;
	_similarityMap[make_pair('z', 'c')] = -3;
	_similarityMap[make_pair('x', 'c')] = -2;
	_similarityMap[make_pair('*', 'c')] = -4;
	_similarityMap[make_pair('a', 'q')] = -1;
	_similarityMap[make_pair('r', 'q')] = 1;
	_similarityMap[make_pair('n', 'q')] = 0;
	_similarityMap[make_pair('d', 'q')] = 0;
	_similarityMap[make_pair('c', 'q')] = -3;
	_similarityMap[make_pair('q', 'q')] = 5;
	_similarityMap[make_pair('e', 'q')] = 2;
	_similarityMap[make_pair('g', 'q')] = -2;
	_similarityMap[make_pair('h', 'q')] = 0;
	_similarityMap[make_pair('i', 'q')] = -3;
	_similarityMap[make_pair('l', 'q')] = -2;
	_similarityMap[make_pair('k', 'q')] = 1;
	_similarityMap[make_pair('m', 'q')] = 0;
	_similarityMap[make_pair('f', 'q')] = -3;
	_similarityMap[make_pair('p', 'q')] = -1;
	_similarityMap[make_pair('s', 'q')] = 0;
	_similarityMap[make_pair('t', 'q')] = -1;
	_similarityMap[make_pair('w', 'q')] = -2;
	_similarityMap[make_pair('y', 'q')] = -1;
	_similarityMap[make_pair('v', 'q')] = -2;
	_similarityMap[make_pair('b', 'q')] = 0;
	_similarityMap[make_pair('z', 'q')] = 3;
	_similarityMap[make_pair('x', 'q')] = -1;
	_similarityMap[make_pair('*', 'q')] = -4;
	_similarityMap[make_pair('a', 'e')] = -1;
	_similarityMap[make_pair('r', 'e')] = 0;
	_similarityMap[make_pair('n', 'e')] = 0;
	_similarityMap[make_pair('d', 'e')] = 2;
	_similarityMap[make_pair('c', 'e')] = -4;
	_similarityMap[make_pair('q', 'e')] = 2;
	_similarityMap[make_pair('e', 'e')] = 5;
	_similarityMap[make_pair('g', 'e')] = -2;
	_similarityMap[make_pair('h', 'e')] = 0;
	_similarityMap[make_pair('i', 'e')] = -3;
	_similarityMap[make_pair('l', 'e')] = -3;
	_similarityMap[make_pair('k', 'e')] = 1;
	_similarityMap[make_pair('m', 'e')] = -2;
	_similarityMap[make_pair('f', 'e')] = -3;
	_similarityMap[make_pair('p', 'e')] = -1;
	_similarityMap[make_pair('s', 'e')] = 0;
	_similarityMap[make_pair('t', 'e')] = -1;
	_similarityMap[make_pair('w', 'e')] = -3;
	_similarityMap[make_pair('y', 'e')] = -2;
	_similarityMap[make_pair('v', 'e')] = -2;
	_similarityMap[make_pair('b', 'e')] = 1;
	_similarityMap[make_pair('z', 'e')] = 4;
	_similarityMap[make_pair('x', 'e')] = -1;
	_similarityMap[make_pair('*', 'e')] = -4;
	_similarityMap[make_pair('a', 'g')] = 0;
	_similarityMap[make_pair('r', 'g')] = -2;
	_similarityMap[make_pair('n', 'g')] = 0;
	_similarityMap[make_pair('d', 'g')] = -1;
	_similarityMap[make_pair('c', 'g')] = -3;
	_similarityMap[make_pair('q', 'g')] = -2;
	_similarityMap[make_pair('e', 'g')] = -2;
	_similarityMap[make_pair('g', 'g')] = 6;
	_similarityMap[make_pair('h', 'g')] = -2;
	_similarityMap[make_pair('i', 'g')] = -4;
	_similarityMap[make_pair('l', 'g')] = -4;
	_similarityMap[make_pair('k', 'g')] = -2;
	_similarityMap[make_pair('m', 'g')] = -3;
	_similarityMap[make_pair('f', 'g')] = -3;
	_similarityMap[make_pair('p', 'g')] = -2;
	_similarityMap[make_pair('s', 'g')] = 0;
	_similarityMap[make_pair('t', 'g')] = -2;
	_similarityMap[make_pair('w', 'g')] = -2;
	_similarityMap[make_pair('y', 'g')] = -3;
	_similarityMap[make_pair('v', 'g')] = -3;
	_similarityMap[make_pair('b', 'g')] = -1;
	_similarityMap[make_pair('z', 'g')] = -2;
	_similarityMap[make_pair('x', 'g')] = -1;
	_similarityMap[make_pair('*', 'g')] = -4;
	_similarityMap[make_pair('a', 'h')] = -2;
	_similarityMap[make_pair('r', 'h')] = 0;
	_similarityMap[make_pair('n', 'h')] = 1;
	_similarityMap[make_pair('d', 'h')] = -1;
	_similarityMap[make_pair('c', 'h')] = -3;
	_similarityMap[make_pair('q', 'h')] = 0;
	_similarityMap[make_pair('e', 'h')] = 0;
	_similarityMap[make_pair('g', 'h')] = -2;
	_similarityMap[make_pair('h', 'h')] = 8;
	_similarityMap[make_pair('i', 'h')] = -3;
	_similarityMap[make_pair('l', 'h')] = -3;
	_similarityMap[make_pair('k', 'h')] = -1;
	_similarityMap[make_pair('m', 'h')] = -2;
	_similarityMap[make_pair('f', 'h')] = -1;
	_similarityMap[make_pair('p', 'h')] = -2;
	_similarityMap[make_pair('s', 'h')] = -1;
	_similarityMap[make_pair('t', 'h')] = -2;
	_similarityMap[make_pair('w', 'h')] = -2;
	_similarityMap[make_pair('y', 'h')] = 2;
	_similarityMap[make_pair('v', 'h')] = -3;
	_similarityMap[make_pair('b', 'h')] = 0;
	_similarityMap[make_pair('z', 'h')] = 0;
	_similarityMap[make_pair('x', 'h')] = -1;
	_similarityMap[make_pair('*', 'h')] = -4;
	_similarityMap[make_pair('a', 'i')] = -1;
	_similarityMap[make_pair('r', 'i')] = -3;
	_similarityMap[make_pair('n', 'i')] = -3;
	_similarityMap[make_pair('d', 'i')] = -3;
	_similarityMap[make_pair('c', 'i')] = -1;
	_similarityMap[make_pair('q', 'i')] = -3;
	_similarityMap[make_pair('e', 'i')] = -3;
	_similarityMap[make_pair('g', 'i')] = -4;
	_similarityMap[make_pair('h', 'i')] = -3;
	_similarityMap[make_pair('i', 'i')] = 4;
	_similarityMap[make_pair('l', 'i')] = 2;
	_similarityMap[make_pair('k', 'i')] = -3;
	_similarityMap[make_pair('m', 'i')] = 1;
	_similarityMap[make_pair('f', 'i')] = 0;
	_similarityMap[make_pair('p', 'i')] = -3;
	_similarityMap[make_pair('s', 'i')] = -2;
	_similarityMap[make_pair('t', 'i')] = -1;
	_similarityMap[make_pair('w', 'i')] = -3;
	_similarityMap[make_pair('y', 'i')] = -1;
	_similarityMap[make_pair('v', 'i')] = 3;
	_similarityMap[make_pair('b', 'i')] = -3;
	_similarityMap[make_pair('z', 'i')] = -3;
	_similarityMap[make_pair('x', 'i')] = -1;
	_similarityMap[make_pair('*', 'i')] = -4;
	_similarityMap[make_pair('a', 'l')] = -1;
	_similarityMap[make_pair('r', 'l')] = -2;
	_similarityMap[make_pair('n', 'l')] = -3;
	_similarityMap[make_pair('d', 'l')] = -4;
	_similarityMap[make_pair('c', 'l')] = -1;
	_similarityMap[make_pair('q', 'l')] = -2;
	_similarityMap[make_pair('e', 'l')] = -3;
	_similarityMap[make_pair('g', 'l')] = -4;
	_similarityMap[make_pair('h', 'l')] = -3;
	_similarityMap[make_pair('i', 'l')] = 2;
	_similarityMap[make_pair('l', 'l')] = 4;
	_similarityMap[make_pair('k', 'l')] = -2;
	_similarityMap[make_pair('m', 'l')] = 2;
	_similarityMap[make_pair('f', 'l')] = 0;
	_similarityMap[make_pair('p', 'l')] = -3;
	_similarityMap[make_pair('s', 'l')] = -2;
	_similarityMap[make_pair('t', 'l')] = -1;
	_similarityMap[make_pair('w', 'l')] = -2;
	_similarityMap[make_pair('y', 'l')] = -1;
	_similarityMap[make_pair('v', 'l')] = 1;
	_similarityMap[make_pair('b', 'l')] = -4;
	_similarityMap[make_pair('z', 'l')] = -3;
	_similarityMap[make_pair('x', 'l')] = -1;
	_similarityMap[make_pair('*', 'l')] = -4;
	_similarityMap[make_pair('a', 'k')] = -1;
	_similarityMap[make_pair('r', 'k')] = 2;
	_similarityMap[make_pair('n', 'k')] = 0;
	_similarityMap[make_pair('d', 'k')] = -1;
	_similarityMap[make_pair('c', 'k')] = -3;
	_similarityMap[make_pair('q', 'k')] = 1;
	_similarityMap[make_pair('e', 'k')] = 1;
	_similarityMap[make_pair('g', 'k')] = -2;
	_similarityMap[make_pair('h', 'k')] = -1;
	_similarityMap[make_pair('i', 'k')] = -3;
	_similarityMap[make_pair('l', 'k')] = -2;
	_similarityMap[make_pair('k', 'k')] = 5;
	_similarityMap[make_pair('m', 'k')] = -1;
	_similarityMap[make_pair('f', 'k')] = -3;
	_similarityMap[make_pair('p', 'k')] = -1;
	_similarityMap[make_pair('s', 'k')] = 0;
	_similarityMap[make_pair('t', 'k')] = -1;
	_similarityMap[make_pair('w', 'k')] = -3;
	_similarityMap[make_pair('y', 'k')] = -2;
	_similarityMap[make_pair('v', 'k')] = -2;
	_similarityMap[make_pair('b', 'k')] = 0;
	_similarityMap[make_pair('z', 'k')] = 1;
	_similarityMap[make_pair('x', 'k')] = -1;
	_similarityMap[make_pair('*', 'k')] = -4;
	_similarityMap[make_pair('a', 'm')] = -1;
	_similarityMap[make_pair('r', 'm')] = -1;
	_similarityMap[make_pair('n', 'm')] = -2;
	_similarityMap[make_pair('d', 'm')] = -3;
	_similarityMap[make_pair('c', 'm')] = -1;
	_similarityMap[make_pair('q', 'm')] = 0;
	_similarityMap[make_pair('e', 'm')] = -2;
	_similarityMap[make_pair('g', 'm')] = -3;
	_similarityMap[make_pair('h', 'm')] = -2;
	_similarityMap[make_pair('i', 'm')] = 1;
	_similarityMap[make_pair('l', 'm')] = 2;
	_similarityMap[make_pair('k', 'm')] = -1;
	_similarityMap[make_pair('m', 'm')] = 5;
	_similarityMap[make_pair('f', 'm')] = 0;
	_similarityMap[make_pair('p', 'm')] = -2;
	_similarityMap[make_pair('s', 'm')] = -1;
	_similarityMap[make_pair('t', 'm')] = -1;
	_similarityMap[make_pair('w', 'm')] = -1;
	_similarityMap[make_pair('y', 'm')] = -1;
	_similarityMap[make_pair('v', 'm')] = 1;
	_similarityMap[make_pair('b', 'm')] = -3;
	_similarityMap[make_pair('z', 'm')] = -1;
	_similarityMap[make_pair('x', 'm')] = -1;
	_similarityMap[make_pair('*', 'm')] = -4;
	_similarityMap[make_pair('a', 'f')] = -2;
	_similarityMap[make_pair('r', 'f')] = -3;
	_similarityMap[make_pair('n', 'f')] = -3;
	_similarityMap[make_pair('d', 'f')] = -3;
	_similarityMap[make_pair('c', 'f')] = -2;
	_similarityMap[make_pair('q', 'f')] = -3;
	_similarityMap[make_pair('e', 'f')] = -3;
	_similarityMap[make_pair('g', 'f')] = -3;
	_similarityMap[make_pair('h', 'f')] = -1;
	_similarityMap[make_pair('i', 'f')] = 0;
	_similarityMap[make_pair('l', 'f')] = 0;
	_similarityMap[make_pair('k', 'f')] = -3;
	_similarityMap[make_pair('m', 'f')] = 0;
	_similarityMap[make_pair('f', 'f')] = 6;
	_similarityMap[make_pair('p', 'f')] = -4;
	_similarityMap[make_pair('s', 'f')] = -2;
	_similarityMap[make_pair('t', 'f')] = -2;
	_similarityMap[make_pair('w', 'f')] = 1;
	_similarityMap[make_pair('y', 'f')] = 3;
	_similarityMap[make_pair('v', 'f')] = -1;
	_similarityMap[make_pair('b', 'f')] = -3;
	_similarityMap[make_pair('z', 'f')] = -3;
	_similarityMap[make_pair('x', 'f')] = -1;
	_similarityMap[make_pair('*', 'f')] = -4;
	_similarityMap[make_pair('a', 'p')] = -1;
	_similarityMap[make_pair('r', 'p')] = -2;
	_similarityMap[make_pair('n', 'p')] = -2;
	_similarityMap[make_pair('d', 'p')] = -1;
	_similarityMap[make_pair('c', 'p')] = -3;
	_similarityMap[make_pair('q', 'p')] = -1;
	_similarityMap[make_pair('e', 'p')] = -1;
	_similarityMap[make_pair('g', 'p')] = -2;
	_similarityMap[make_pair('h', 'p')] = -2;
	_similarityMap[make_pair('i', 'p')] = -3;
	_similarityMap[make_pair('l', 'p')] = -3;
	_similarityMap[make_pair('k', 'p')] = -1;
	_similarityMap[make_pair('m', 'p')] = -2;
	_similarityMap[make_pair('f', 'p')] = -4;
	_similarityMap[make_pair('p', 'p')] = 7;
	_similarityMap[make_pair('s', 'p')] = -1;
	_similarityMap[make_pair('t', 'p')] = -1;
	_similarityMap[make_pair('w', 'p')] = -4;
	_similarityMap[make_pair('y', 'p')] = -3;
	_similarityMap[make_pair('v', 'p')] = -2;
	_similarityMap[make_pair('b', 'p')] = -2;
	_similarityMap[make_pair('z', 'p')] = -1;
	_similarityMap[make_pair('x', 'p')] = -2;
	_similarityMap[make_pair('*', 'p')] = -4;
	_similarityMap[make_pair('a', 's')] = 1;
	_similarityMap[make_pair('r', 's')] = -1;
	_similarityMap[make_pair('n', 's')] = 1;
	_similarityMap[make_pair('d', 's')] = 0;
	_similarityMap[make_pair('c', 's')] = -1;
	_similarityMap[make_pair('q', 's')] = 0;
	_similarityMap[make_pair('e', 's')] = 0;
	_similarityMap[make_pair('g', 's')] = 0;
	_similarityMap[make_pair('h', 's')] = -1;
	_similarityMap[make_pair('i', 's')] = -2;
	_similarityMap[make_pair('l', 's')] = -2;
	_similarityMap[make_pair('k', 's')] = 0;
	_similarityMap[make_pair('m', 's')] = -1;
	_similarityMap[make_pair('f', 's')] = -2;
	_similarityMap[make_pair('p', 's')] = -1;
	_similarityMap[make_pair('s', 's')] = 4;
	_similarityMap[make_pair('t', 's')] = 1;
	_similarityMap[make_pair('w', 's')] = -3;
	_similarityMap[make_pair('y', 's')] = -2;
	_similarityMap[make_pair('v', 's')] = -2;
	_similarityMap[make_pair('b', 's')] = 0;
	_similarityMap[make_pair('z', 's')] = 0;
	_similarityMap[make_pair('x', 's')] = 0;
	_similarityMap[make_pair('*', 's')] = -4;
	_similarityMap[make_pair('a', 't')] = 0;
	_similarityMap[make_pair('r', 't')] = -1;
	_similarityMap[make_pair('n', 't')] = 0;
	_similarityMap[make_pair('d', 't')] = -1;
	_similarityMap[make_pair('c', 't')] = -1;
	_similarityMap[make_pair('q', 't')] = -1;
	_similarityMap[make_pair('e', 't')] = -1;
	_similarityMap[make_pair('g', 't')] = -2;
	_similarityMap[make_pair('h', 't')] = -2;
	_similarityMap[make_pair('i', 't')] = -1;
	_similarityMap[make_pair('l', 't')] = -1;
	_similarityMap[make_pair('k', 't')] = -1;
	_similarityMap[make_pair('m', 't')] = -1;
	_similarityMap[make_pair('f', 't')] = -2;
	_similarityMap[make_pair('p', 't')] = -1;
	_similarityMap[make_pair('s', 't')] = 1;
	_similarityMap[make_pair('t', 't')] = 5;
	_similarityMap[make_pair('w', 't')] = -2;
	_similarityMap[make_pair('y', 't')] = -2;
	_similarityMap[make_pair('v', 't')] = 0;
	_similarityMap[make_pair('b', 't')] = -1;
	_similarityMap[make_pair('z', 't')] = -1;
	_similarityMap[make_pair('x', 't')] = 0;
	_similarityMap[make_pair('*', 't')] = -4;
	_similarityMap[make_pair('a', 'w')] = -3;
	_similarityMap[make_pair('r', 'w')] = -3;
	_similarityMap[make_pair('n', 'w')] = -4;
	_similarityMap[make_pair('d', 'w')] = -4;
	_similarityMap[make_pair('c', 'w')] = -2;
	_similarityMap[make_pair('q', 'w')] = -2;
	_similarityMap[make_pair('e', 'w')] = -3;
	_similarityMap[make_pair('g', 'w')] = -2;
	_similarityMap[make_pair('h', 'w')] = -2;
	_similarityMap[make_pair('i', 'w')] = -3;
	_similarityMap[make_pair('l', 'w')] = -2;
	_similarityMap[make_pair('k', 'w')] = -3;
	_similarityMap[make_pair('m', 'w')] = -1;
	_similarityMap[make_pair('f', 'w')] = 1;
	_similarityMap[make_pair('p', 'w')] = -4;
	_similarityMap[make_pair('s', 'w')] = -3;
	_similarityMap[make_pair('t', 'w')] = -2;
	_similarityMap[make_pair('w', 'w')] = 11;
	_similarityMap[make_pair('y', 'w')] = 2;
	_similarityMap[make_pair('v', 'w')] = -3;
	_similarityMap[make_pair('b', 'w')] = -4;
	_similarityMap[make_pair('z', 'w')] = -3;
	_similarityMap[make_pair('x', 'w')] = -2;
	_similarityMap[make_pair('*', 'w')] = -4;
	_similarityMap[make_pair('a', 'y')] = -2;
	_similarityMap[make_pair('r', 'y')] = -2;
	_similarityMap[make_pair('n', 'y')] = -2;
	_similarityMap[make_pair('d', 'y')] = -3;
	_similarityMap[make_pair('c', 'y')] = -2;
	_similarityMap[make_pair('q', 'y')] = -1;
	_similarityMap[make_pair('e', 'y')] = -2;
	_similarityMap[make_pair('g', 'y')] = -3;
	_similarityMap[make_pair('h', 'y')] = 2;
	_similarityMap[make_pair('i', 'y')] = -1;
	_similarityMap[make_pair('l', 'y')] = -1;
	_similarityMap[make_pair('k', 'y')] = -2;
	_similarityMap[make_pair('m', 'y')] = -1;
	_similarityMap[make_pair('f', 'y')] = 3;
	_similarityMap[make_pair('p', 'y')] = -3;
	_similarityMap[make_pair('s', 'y')] = -2;
	_similarityMap[make_pair('t', 'y')] = -2;
	_similarityMap[make_pair('w', 'y')] = 2;
	_similarityMap[make_pair('y', 'y')] = 7;
	_similarityMap[make_pair('v', 'y')] = -1;
	_similarityMap[make_pair('b', 'y')] = -3;
	_similarityMap[make_pair('z', 'y')] = -2;
	_similarityMap[make_pair('x', 'y')] = -1;
	_similarityMap[make_pair('*', 'y')] = -4;
	_similarityMap[make_pair('a', 'v')] = 0;
	_similarityMap[make_pair('r', 'v')] = -3;
	_similarityMap[make_pair('n', 'v')] = -3;
	_similarityMap[make_pair('d', 'v')] = -3;
	_similarityMap[make_pair('c', 'v')] = -1;
	_similarityMap[make_pair('q', 'v')] = -2;
	_similarityMap[make_pair('e', 'v')] = -2;
	_similarityMap[make_pair('g', 'v')] = -3;
	_similarityMap[make_pair('h', 'v')] = -3;
	_similarityMap[make_pair('i', 'v')] = 3;
	_similarityMap[make_pair('l', 'v')] = 1;
	_similarityMap[make_pair('k', 'v')] = -2;
	_similarityMap[make_pair('m', 'v')] = 1;
	_similarityMap[make_pair('f', 'v')] = -1;
	_similarityMap[make_pair('p', 'v')] = -2;
	_similarityMap[make_pair('s', 'v')] = -2;
	_similarityMap[make_pair('t', 'v')] = 0;
	_similarityMap[make_pair('w', 'v')] = -3;
	_similarityMap[make_pair('y', 'v')] = -1;
	_similarityMap[make_pair('v', 'v')] = 4;
	_similarityMap[make_pair('b', 'v')] = -3;
	_similarityMap[make_pair('z', 'v')] = -2;
	_similarityMap[make_pair('x', 'v')] = -1;
	_similarityMap[make_pair('*', 'v')] = -4;
	_similarityMap[make_pair('a', 'b')] = -2;
	_similarityMap[make_pair('r', 'b')] = -1;
	_similarityMap[make_pair('n', 'b')] = 3;
	_similarityMap[make_pair('d', 'b')] = 4;
	_similarityMap[make_pair('c', 'b')] = -3;
	_similarityMap[make_pair('q', 'b')] = 0;
	_similarityMap[make_pair('e', 'b')] = 1;
	_similarityMap[make_pair('g', 'b')] = -1;
	_similarityMap[make_pair('h', 'b')] = 0;
	_similarityMap[make_pair('i', 'b')] = -3;
	_similarityMap[make_pair('l', 'b')] = -4;
	_similarityMap[make_pair('k', 'b')] = 0;
	_similarityMap[make_pair('m', 'b')] = -3;
	_similarityMap[make_pair('f', 'b')] = -3;
	_similarityMap[make_pair('p', 'b')] = -2;
	_similarityMap[make_pair('s', 'b')] = 0;
	_similarityMap[make_pair('t', 'b')] = -1;
	_similarityMap[make_pair('w', 'b')] = -4;
	_similarityMap[make_pair('y', 'b')] = -3;
	_similarityMap[make_pair('v', 'b')] = -3;
	_similarityMap[make_pair('b', 'b')] = 4;
	_similarityMap[make_pair('z', 'b')] = 1;
	_similarityMap[make_pair('x', 'b')] = -1;
	_similarityMap[make_pair('*', 'b')] = -4;
	_similarityMap[make_pair('a', 'z')] = -1;
	_similarityMap[make_pair('r', 'z')] = 0;
	_similarityMap[make_pair('n', 'z')] = 0;
	_similarityMap[make_pair('d', 'z')] = 1;
	_similarityMap[make_pair('c', 'z')] = -3;
	_similarityMap[make_pair('q', 'z')] = 3;
	_similarityMap[make_pair('e', 'z')] = 4;
	_similarityMap[make_pair('g', 'z')] = -2;
	_similarityMap[make_pair('h', 'z')] = 0;
	_similarityMap[make_pair('i', 'z')] = -3;
	_similarityMap[make_pair('l', 'z')] = -3;
	_similarityMap[make_pair('k', 'z')] = 1;
	_similarityMap[make_pair('m', 'z')] = -1;
	_similarityMap[make_pair('f', 'z')] = -3;
	_similarityMap[make_pair('p', 'z')] = -1;
	_similarityMap[make_pair('s', 'z')] = 0;
	_similarityMap[make_pair('t', 'z')] = -1;
	_similarityMap[make_pair('w', 'z')] = -3;
	_similarityMap[make_pair('y', 'z')] = -2;
	_similarityMap[make_pair('v', 'z')] = -2;
	_similarityMap[make_pair('b', 'z')] = 1;
	_similarityMap[make_pair('z', 'z')] = 4;
	_similarityMap[make_pair('x', 'z')] = -1;
	_similarityMap[make_pair('*', 'z')] = -4;
	_similarityMap[make_pair('a', 'x')] = 0;
	_similarityMap[make_pair('r', 'x')] = -1;
	_similarityMap[make_pair('n', 'x')] = -1;
	_similarityMap[make_pair('d', 'x')] = -1;
	_similarityMap[make_pair('c', 'x')] = -2;
	_similarityMap[make_pair('q', 'x')] = -1;
	_similarityMap[make_pair('e', 'x')] = -1;
	_similarityMap[make_pair('g', 'x')] = -1;
	_similarityMap[make_pair('h', 'x')] = -1;
	_similarityMap[make_pair('i', 'x')] = -1;
	_similarityMap[make_pair('l', 'x')] = -1;
	_similarityMap[make_pair('k', 'x')] = -1;
	_similarityMap[make_pair('m', 'x')] = -1;
	_similarityMap[make_pair('f', 'x')] = -1;
	_similarityMap[make_pair('p', 'x')] = -2;
	_similarityMap[make_pair('s', 'x')] = 0;
	_similarityMap[make_pair('t', 'x')] = 0;
	_similarityMap[make_pair('w', 'x')] = -2;
	_similarityMap[make_pair('y', 'x')] = -1;
	_similarityMap[make_pair('v', 'x')] = -1;
	_similarityMap[make_pair('b', 'x')] = -1;
	_similarityMap[make_pair('z', 'x')] = -1;
	_similarityMap[make_pair('x', 'x')] = -1;
	_similarityMap[make_pair('*', 'x')] = -4;
	_similarityMap[make_pair('a', '*')] = -4;
	_similarityMap[make_pair('r', '*')] = -4;
	_similarityMap[make_pair('n', '*')] = -4;
	_similarityMap[make_pair('d', '*')] = -4;
	_similarityMap[make_pair('c', '*')] = -4;
	_similarityMap[make_pair('q', '*')] = -4;
	_similarityMap[make_pair('e', '*')] = -4;
	_similarityMap[make_pair('g', '*')] = -4;
	_similarityMap[make_pair('h', '*')] = -4;
	_similarityMap[make_pair('i', '*')] = -4;
	_similarityMap[make_pair('l', '*')] = -4;
	_similarityMap[make_pair('k', '*')] = -4;
	_similarityMap[make_pair('m', '*')] = -4;
	_similarityMap[make_pair('f', '*')] = -4;
	_similarityMap[make_pair('p', '*')] = -4;
	_similarityMap[make_pair('s', '*')] = -4;
	_similarityMap[make_pair('t', '*')] = -4;
	_similarityMap[make_pair('w', '*')] = -4;
	_similarityMap[make_pair('y', '*')] = -4;
	_similarityMap[make_pair('v', '*')] = -4;
	_similarityMap[make_pair('b', '*')] = -4;
	_similarityMap[make_pair('z', '*')] = -4;
	_similarityMap[make_pair('x', '*')] = -4;
	_similarityMap[make_pair('*', '*')] = 1;
}

int read_input_seqs::get_pair_score(char a, char b)
{
	char lower_a = tolower(a);
	char lower_b = tolower(b);
	pair<char, char> curr_pair = make_pair(lower_a, lower_b);

	if (_similarityMap.find(curr_pair) == _similarityMap.end())
	{
		// not in map
		if (a == b)
		{
			return _default_match_score;
		}
		else
		{
			return _default_mismatch_score;
		}
	}
	// exists in map
	int value = _similarityMap[curr_pair];
	return value;

}