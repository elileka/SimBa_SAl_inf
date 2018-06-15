// code by Eli Levy Karin

#include "read_chop_tables.h"

double read_chop_tables::get_chop_prob(const char chop_type, size_t i, size_t j)
{
	pair<size_t, size_t> curr_ij = make_pair(i, j);
	size_t num_to_divide_by = (i + j) > 680 ? 680 : (i + j); // to avoid too small numbers...
	double value_to_return_if_not_found = _min_probability / 100.0; // that is very very small but not 0
	if (chop_type == 'N')
	{
		if (_N_chop_probs_map.find(curr_ij) == _N_chop_probs_map.end())
		{
			// not in map
			return(value_to_return_if_not_found / exp(num_to_divide_by)); // smaller if i and j are big
		}
		else
		{
			// exists in map
			 return(_N_chop_probs_map[curr_ij]);
		}
	}
	else if (chop_type == 'L')
	{
		if (_L_chop_probs_map.find(curr_ij) == _L_chop_probs_map.end())
		{
			// not in map
			return(value_to_return_if_not_found / exp(num_to_divide_by + 5)); // smaller if i and j are big, L is less probable than N
		}
		else
		{
			// exists in map
			return(_L_chop_probs_map[curr_ij]);
		}
	}
	else if (chop_type == 'R')
	{
		if (_R_chop_probs_map.find(curr_ij) == _R_chop_probs_map.end())
		{
			// not in map
			return(value_to_return_if_not_found / exp(num_to_divide_by + 5)); // smaller if i and j are big, R is less probable than N
		}
		else
		{
			// exists in map
			return(_R_chop_probs_map[curr_ij]);
		}
	}
	else if (chop_type == 'B')
	{
		if (_B_chop_probs_map.find(curr_ij) == _B_chop_probs_map.end())
		{
			// not in map
			return(value_to_return_if_not_found / exp(num_to_divide_by + 10)); // smaller if i and j are big, B is less probable than N,L,R
		}
		else
		{
			// exists in map
			return(_B_chop_probs_map[curr_ij]);
		}
	}
	else
	{
		cerr << "illegal chop type: " << chop_type << ", allowed types: [N,L,R,B]" << endl;
		exit(1);
	}
	
	// we do not get here:
	return -1.0;
}

void read_chop_tables::read_tables_from_file()
{
	ifstream maps_file(_chop_tables_file_string);
	string line;

	getline(maps_file, line); // get rid of header line
	size_t num_chop_lines_read = 0;
	while (getline(maps_file, line))
	{
		stringstream linestream(line);
		// type_chop	i	j	num_times_chop_observed	total_num_chops_of_this_type	chop_probability
		string chop_type; // N, L, R, B
		size_t i;
		size_t j;
		size_t num_times_chop_observed;
		size_t total_num_chops_of_this_type;
		double chop_prob;

		getline(linestream, chop_type, '\t'); // read up-to the first tab (discard tab)
		linestream >> i >> j >> num_times_chop_observed >> total_num_chops_of_this_type >> chop_prob; // read the size_t using the operator >>
		
		if (_min_probability > chop_prob)
		{
			_min_probability = chop_prob;
		}

		pair<size_t, size_t> curr_ij = make_pair(i, j);
		if (chop_type.compare("N") == 0)
		{
			_N_chop_probs_map[curr_ij] = chop_prob;
		}
		else if (chop_type.compare("L") == 0)
		{
			_L_chop_probs_map[curr_ij] = chop_prob;
			if (i > _max_i_in_L_table)
			{
				_max_i_in_L_table = i;
			}
			if (j > _max_j_in_L_table)
			{
				_max_j_in_L_table = j;
			}
		}
		else if (chop_type.compare("R") == 0)
		{
			_R_chop_probs_map[curr_ij] = chop_prob;
			if (i > _max_i_in_R_table)
			{
				_max_i_in_R_table = i;
			}
			if (j > _max_j_in_R_table)
			{
				_max_j_in_R_table = j;
			}
		}
		else if (chop_type.compare("B") == 0)
		{
			_B_chop_probs_map[curr_ij] = chop_prob;
		}
		else
		{
			cerr << "unrecognized type of chop! " << chop_type << endl;
			exit(1);
		}
		num_chop_lines_read++;
	}
	if (num_chop_lines_read == 0)
	{
		cerr << "something wrong with reading the table - no chop lines read from: " << _chop_tables_file_string << endl;
		exit(1);
	}
}
