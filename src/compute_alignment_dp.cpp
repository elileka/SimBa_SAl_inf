// code by Eli Levy Karin

#include "compute_alignment_dp.h"

vector<string> compute_alignment_dp::get_sampled_alignment(double & sampled_alignment_log_prob, bool should_take_max, int band_width)
{
	// fill the S table with the sampled path info
	fill_S_table_corner_cutting(should_take_max, band_width);
	chop_prob sampled_last_chop = get_sampled_path_last_chop(should_take_max, sampled_alignment_log_prob);

	vector<size_t> anc_seq = _coded_seqs[0];
	vector<size_t> des_seq = _coded_seqs[1];

	vector<vector<size_t>> sampled_alignment;
	vector<string> sampled_alignment_strings;

	if (sampled_last_chop.get_chop_type() == 'B') // non-homologous alignment
	{
		for (size_t n = 0; n < _length_of_anc; n++)
		{
			vector<size_t> curr_pos;
			curr_pos.push_back(anc_seq[n]);
			curr_pos.push_back(0); // anc deleted from des
			sampled_alignment.push_back(curr_pos);
		}
		for (size_t m = 0; m < _length_of_des; m++)
		{
			vector<size_t> curr_pos;
			curr_pos.push_back(0); // inserted in des
			curr_pos.push_back(des_seq[m]);
			sampled_alignment.push_back(curr_pos);
		}

		//string sampled_alignment_anc_str = get_str_aligned_seq(0, sampled_alignment);
		//string sampled_alignment_des_str = get_str_aligned_seq(1, sampled_alignment);
		string sampled_alignment_anc_str = _read_input_seqs_obj.get_str_from_coded_alignment(0, sampled_alignment);
		string sampled_alignment_des_str = _read_input_seqs_obj.get_str_from_coded_alignment(1, sampled_alignment);

		sampled_alignment_strings.push_back(sampled_alignment_anc_str);
		sampled_alignment_strings.push_back(sampled_alignment_des_str);
		return sampled_alignment_strings;
	}

	// if we are here, the last chop is of type R
	// sampled_last_chop.print_chop_info(); // for debug

	vector<vector<size_t>> sampled_alignment_reversed;
	int anc_ind = _length_of_anc - 1; // initialize
	int des_ind = _length_of_des - 1; // initialize
	for (size_t m = 0; m < sampled_last_chop.get_j(); m++)
	{
		vector<size_t> curr_pos;
		curr_pos.push_back(0);
		curr_pos.push_back(des_seq[des_ind]); // inserted in des
		sampled_alignment_reversed.push_back(curr_pos);
		des_ind--;
	}
	for (size_t n = 0; n < sampled_last_chop.get_i(); n++)
	{
		vector<size_t> curr_pos;
		curr_pos.push_back(anc_seq[anc_ind]);
		curr_pos.push_back(0); // anc deleted from des
		sampled_alignment_reversed.push_back(curr_pos);
		anc_ind--;
	}

	while ((anc_ind >= 0) && (des_ind >= 0))
	{
		chop_prob curr_chop = _S_table[anc_ind][des_ind].first; // either L or N type!
																// curr_chop.print_chop_info(); // for debug

																// push the matched position:
		vector<size_t> curr_pos;
		curr_pos.push_back(anc_seq[anc_ind]);
		curr_pos.push_back(des_seq[des_ind]);
		sampled_alignment_reversed.push_back(curr_pos);
		anc_ind--;
		des_ind--;

		// push the non-homologous positions
		for (size_t m = 0; m < curr_chop.get_j(); m++)
		{
			vector<size_t> curr_pos;
			curr_pos.push_back(0);
			curr_pos.push_back(des_seq[des_ind]); // inserted in des
			sampled_alignment_reversed.push_back(curr_pos);
			des_ind--;
		}
		for (size_t n = 0; n < curr_chop.get_i(); n++)
		{
			vector<size_t> curr_pos;
			curr_pos.push_back(anc_seq[anc_ind]);
			curr_pos.push_back(0); // anc deleted from des
			sampled_alignment_reversed.push_back(curr_pos);
			anc_ind--;
		}
	}

	// if some ancestral chars are left:
	while (anc_ind >= 0)
	{
		vector<size_t> curr_pos;
		curr_pos.push_back(anc_seq[anc_ind]);
		curr_pos.push_back(0); // anc deleted from des
		sampled_alignment_reversed.push_back(curr_pos);
		anc_ind--;
	}

	// if some descendant chars are left:
	while (des_ind >= 0)
	{
		vector<size_t> curr_pos;
		curr_pos.push_back(0);
		curr_pos.push_back(des_seq[des_ind]); // inserted in des
		sampled_alignment_reversed.push_back(curr_pos);
		des_ind--;
	}

	// reverse the order
	size_t length_of_alignment = sampled_alignment_reversed.size();
	for (size_t alignment_ind = 0; alignment_ind < length_of_alignment; alignment_ind++)
	{
		sampled_alignment.push_back(sampled_alignment_reversed[length_of_alignment - 1 - alignment_ind]);
	}


	//string sampled_alignment_anc_str = get_str_aligned_seq(0, sampled_alignment);
	//string sampled_alignment_des_str = get_str_aligned_seq(1, sampled_alignment);
	string sampled_alignment_anc_str = _read_input_seqs_obj.get_str_from_coded_alignment(0, sampled_alignment);
	string sampled_alignment_des_str = _read_input_seqs_obj.get_str_from_coded_alignment(1, sampled_alignment);

	sampled_alignment_strings.push_back(sampled_alignment_anc_str);
	sampled_alignment_strings.push_back(sampled_alignment_des_str);
	return sampled_alignment_strings;
}


size_t compute_alignment_dp::sample_index_gumbel_max(vector<double>& log_probs)
{
	// this methd implements the Gumble-max sampling to allow
	// sampling from a vector of log-probabilities while
	// working in log-space
	// references:
	// https://stats.stackexchange.com/questions/64081/how-do-i-sample-from-a-discrete-categorical-distribution-in-log-space
	// https://en.wikipedia.org/wiki/Categorical_distribution#Sampling_via_the_Gumbel_distribution

	std::random_device rd;  // will be used to obtain a seed for the random number engine
	std::mt19937 gen(rd()); // standard mersenne_twister_engine seeded with rd()
	std::uniform_real_distribution<> dis(0.0, 1.0); // random between 0 and 1

	size_t chosen_index = 0;
	double max_gi_plus_log_prob_i = 0;

	for (size_t vec_ind = 0; vec_ind < log_probs.size(); vec_ind++)
	{
		double uniform_helper = dis(gen); // this will give us the required randomness
		double curr_gi = -log(-log(uniform_helper));

		if (vec_ind == 0) // first iteration
		{
			max_gi_plus_log_prob_i = curr_gi + log_probs[vec_ind];
			chosen_index = 0;
		}
		else if (max_gi_plus_log_prob_i < (curr_gi + log_probs[vec_ind]))
		{
			max_gi_plus_log_prob_i = curr_gi + log_probs[vec_ind];
			chosen_index = vec_ind;
		}
	}

	return (chosen_index);
}

void compute_alignment_dp::get_combs(int band_width, vector<vector<size_t>>& pairs_to_consider)
{
	/*
	pairs_to_consider is a data structure to hold <i,j> combinations that are
	worth considering.
	The iteration on this vector assures the right order of looking at the combinations.

	if band_width == -1 --> all combinations will be considered (no corner cutting)
	*/
	cout << "going to use band_width = (" << band_width << ") to construct the pairs to consider" << endl; // for debug
	cout << "_length_of_des is: " << _length_of_des << endl; // for debug
	cout << "_length_of_anc is: " << _length_of_anc << endl; // for debug

	size_t num_pairs_to_consider = 0; // for debug

	vector<vector<vector<size_t>>> i_and_j_combs_and_prev_combs;
	for (size_t anc_ind_i = 0; anc_ind_i < _length_of_anc; anc_ind_i++)
	{
		size_t des_min_ind = 0; // border without band
		size_t des_max_ind = _length_of_des - 1; // border without band
												 // deal with band borders if there is a band:
		if (band_width > -1)
		{
			size_t diff_of_lengths = 0;
			if (_length_of_anc > _length_of_des)
			{
				diff_of_lengths = _length_of_anc - _length_of_des;

				if (anc_ind_i > diff_of_lengths + (size_t)band_width)
				{
					des_min_ind = anc_ind_i - diff_of_lengths - (size_t)band_width;
				}
				if (des_max_ind > anc_ind_i + (size_t)band_width)
				{
					des_max_ind = anc_ind_i + (size_t)band_width;
				}
			}
			else
			{
				diff_of_lengths = _length_of_des - _length_of_anc;
				if (anc_ind_i > (size_t)band_width)
				{
					des_min_ind = anc_ind_i - (size_t)band_width;
				}
				if (des_max_ind > anc_ind_i + diff_of_lengths + (size_t)band_width)
				{
					des_max_ind = anc_ind_i + diff_of_lengths + (size_t)band_width;
				}
			}
		}
		// end deal with band borders

		// cout << "anc_ind_i: " << anc_ind_i << " with: des_ind_j = [" << des_min_ind << "," << des_max_ind << "]" << endl;

		for (size_t des_ind_j = des_min_ind; des_ind_j <= des_max_ind; des_ind_j++)
		{
			vector<vector<size_t>> curr_comb_and_its_prevs;

			vector<size_t> curr_combinaion;
			curr_combinaion.push_back(anc_ind_i);
			curr_combinaion.push_back(des_ind_j);
			num_pairs_to_consider++; // for debug

			pairs_to_consider.push_back(curr_combinaion);
		}
	}

	cout << "in total: " << num_pairs_to_consider << " pairs will be considered" << endl;
}

void compute_alignment_dp::get_prev_combs_for_pair(int band_width, vector<vector<size_t>>& prev_pairs_to_consider, size_t anc_ind_i, size_t des_ind_j)
{
	for (size_t prev_anc_ind_i = 0; prev_anc_ind_i < anc_ind_i; prev_anc_ind_i++)
	{
		size_t prev_des_min_ind = 0; // border without band
		size_t prev_des_max_ind_plus_1 = des_ind_j; // border without band
													// deal with band borders if there is a band:
		if (band_width > -1)
		{
			size_t diff_of_lengths = 0;
			if (_length_of_anc > _length_of_des)
			{
				diff_of_lengths = _length_of_anc - _length_of_des;
				if (prev_anc_ind_i > diff_of_lengths + (size_t)band_width)
				{
					prev_des_min_ind = prev_anc_ind_i - diff_of_lengths - (size_t)band_width;
				}
				if (prev_des_max_ind_plus_1 > prev_anc_ind_i + (size_t)band_width + 1)
				{
					prev_des_max_ind_plus_1 = prev_anc_ind_i + (size_t)band_width + 1;
				}
			}
			else
			{
				diff_of_lengths = _length_of_des - _length_of_anc;
				if (prev_anc_ind_i > (size_t)band_width)
				{
					prev_des_min_ind = prev_anc_ind_i - (size_t)band_width;
				}
				if (prev_des_max_ind_plus_1 > prev_anc_ind_i + diff_of_lengths + (size_t)band_width + 1)
				{
					prev_des_max_ind_plus_1 = prev_anc_ind_i + diff_of_lengths + (size_t)band_width + 1;
				}
			}
		}
		// end deal with band borders

		for (size_t prev_des_ind_j = prev_des_min_ind; prev_des_ind_j < prev_des_max_ind_plus_1; prev_des_ind_j++)
		{
			vector<size_t> prev_combinaion;
			prev_combinaion.push_back(prev_anc_ind_i);
			prev_combinaion.push_back(prev_des_ind_j);
			prev_pairs_to_consider.push_back(prev_combinaion);
		}
	}
}

void compute_alignment_dp::get_backward_combs(int band_width, vector<vector<size_t>>& pairs_to_consider_backward)
{
	/*
	This function relies on the technique of get_combs
	by mapping the original indeces to reversed ones :-)
	*/

	vector<vector<size_t>> i_and_j_combs;
	get_combs(band_width, i_and_j_combs);

	for (size_t i_and_j_comb_ind = 0; i_and_j_comb_ind < i_and_j_combs.size(); i_and_j_comb_ind++)
	{
		size_t anc_ind_i = i_and_j_combs[i_and_j_comb_ind][0];
		size_t des_ind_j = i_and_j_combs[i_and_j_comb_ind][1];

		size_t anc_ind_i_backward = _length_of_anc - anc_ind_i - 1; // backward translation
		size_t des_ind_j_backward = _length_of_des - des_ind_j - 1; // backward translation

		vector<size_t> curr_combinaion_backward;
		curr_combinaion_backward.push_back(anc_ind_i_backward);
		curr_combinaion_backward.push_back(des_ind_j_backward);
		pairs_to_consider_backward.push_back(curr_combinaion_backward);
	}
}

void compute_alignment_dp::get_backward_prev_combs_for_pair(int band_width, vector<vector<size_t>>& prev_pairs_to_consider_backward, size_t anc_ind_i, size_t des_ind_j)
{
	/*
	This function relies on the technique of get_prev_combs_for_pair
	by mapping the original indeces to reversed ones :-)
	*/

	size_t forward_anc_ind_i = _length_of_anc - anc_ind_i - 1;
	size_t forward_des_ind_j = _length_of_des - des_ind_j - 1;

	vector<vector<size_t>> prev_i_and_j_combs;
	get_prev_combs_for_pair(band_width, prev_i_and_j_combs, forward_anc_ind_i, forward_des_ind_j);

	for (size_t prev_i_j_combs_ind = 1; prev_i_j_combs_ind < prev_i_and_j_combs.size(); prev_i_j_combs_ind++)
	{
		size_t curr_anc_prev_ind = prev_i_and_j_combs[prev_i_j_combs_ind][0];
		size_t curr_des_prev_ind = prev_i_and_j_combs[prev_i_j_combs_ind][1];

		size_t curr_anc_prev_ind_backward = _length_of_anc - curr_anc_prev_ind - 1; // backward translation
		size_t curr_des_prev_ind_backward = _length_of_des - curr_des_prev_ind - 1; // backward translation

		vector<size_t> curr_prev_combinaion_backward;
		curr_prev_combinaion_backward.push_back(curr_anc_prev_ind_backward);
		curr_prev_combinaion_backward.push_back(curr_des_prev_ind_backward);
		prev_pairs_to_consider_backward.push_back(curr_prev_combinaion_backward);
	}
}

chop_prob compute_alignment_dp::get_sampled_path_last_chop(bool should_take_max, double & sampled_alignment_log_prob)
{
	vector<size_t> anc_seq = _coded_seqs[0];
	vector<size_t> des_seq = _coded_seqs[1];

	chop_prob the_B_chop(_chop_tables_obj, 'B', _length_of_anc, _length_of_des, 0, 0, des_seq, _branch_length_t, _is_jc, _quick_jtt_obj);

	vector<chop_prob> last_chop_options; // can be used for sampling, if needed
	vector<double> alignments_log_probs; // can be used for sampling, if needed
	last_chop_options.push_back(the_B_chop); // can be used for sampling, if needed
	alignments_log_probs.push_back(the_B_chop.get_chop_log_prob()); // can be used for sampling, if needed

	double best_alignment_log_prob = the_B_chop.get_chop_log_prob(); // maximum alignment, if needed
	chop_prob best_last_chop = the_B_chop;

	for (size_t n = 0; n < _length_of_anc; n++)
	{
		for (size_t m = 0; m < _length_of_des; m++)
		{
			vector<size_t> inserted_des_chars_internal; // in case m = 0, this will be empty
			if (m > 0)
			{
				for (size_t k = (_length_of_des - m); k < _length_of_des; k++)
				{
					inserted_des_chars_internal.push_back(des_seq[k]);
				}
			}
			chop_prob curr_Rnm_chop(_chop_tables_obj, 'R', n, m, 0, 0, inserted_des_chars_internal, _branch_length_t, _is_jc, _quick_jtt_obj);
			double curr_Rnm_total_log_prob = curr_Rnm_chop.get_chop_log_prob();
			double curr_prev_S_prop = _S_table[_length_of_anc - n - 1][_length_of_des - m - 1].second;

			if (curr_prev_S_prop > 0) // this cell was not computed due to corner cutting
			{
				continue; // we don't consider this combination!
			}

			double curr_alignment_log_prob = (curr_prev_S_prop + curr_Rnm_total_log_prob);

			last_chop_options.push_back(curr_Rnm_chop); // can be used for sampling, if needed
			alignments_log_probs.push_back(curr_alignment_log_prob); // can be used for sampling, if needed

			if (best_alignment_log_prob < curr_alignment_log_prob)
			{
				best_alignment_log_prob = curr_alignment_log_prob;
				best_last_chop = curr_Rnm_chop;
			}
		}
	}

	if (should_take_max) // in case the maximum alignment is desired
	{
		sampled_alignment_log_prob = best_alignment_log_prob;
		return (best_last_chop);
	}
	else // sampling is required
	{
		size_t sampled_ind = sample_index_gumbel_max(alignments_log_probs);
		sampled_alignment_log_prob = alignments_log_probs[sampled_ind];
		return (last_chop_options[sampled_ind]);
	}

}

double compute_alignment_dp::compute_conditional_alignment_total_log_prob(int band_width)
{
	// fill in _P_table - if this is the first call:
	if (!_is_P_computed)
	{
		fill_P_table_corner_cutting(band_width);
	}

	// now compute according to eq 15 of Miklos et al 2004 with zero indices:

	vector<size_t> anc_seq = _coded_seqs[0];
	vector<size_t> des_seq = _coded_seqs[1];

	chop_prob the_B_chop(_chop_tables_obj, 'B', _length_of_anc, _length_of_des, 0, 0, des_seq, _branch_length_t, _is_jc, _quick_jtt_obj);
	double curr_B_total_log_prob = the_B_chop.get_chop_log_prob();

	vector<double> log_probs; // will be used for summing probs
	log_probs.push_back(curr_B_total_log_prob); // add the first log prob
	double max_obs_log_prob = curr_B_total_log_prob; // initialize

	for (size_t n = 0; n < _length_of_anc; n++)
	{
		for (size_t m = 0; m < _length_of_des; m++)
		{
			vector<size_t> inserted_des_chars_internal; // in case m = 0, this will be empty
			if (m > 0)
			{
				for (size_t k = (_length_of_des - m); k < _length_of_des; k++)
				{
					inserted_des_chars_internal.push_back(des_seq[k]);
				}
			}
			chop_prob curr_Rnm_chop(_chop_tables_obj, 'R', n, m, 0, 0, inserted_des_chars_internal, _branch_length_t, _is_jc, _quick_jtt_obj);
			double curr_Rnm_total_log_prob = curr_Rnm_chop.get_chop_log_prob();
			double curr_prev_P_log = _P_table[_length_of_anc - n - 1][_length_of_des - m - 1];

			if (curr_prev_P_log > 0) // illegal value - this means cell is irrelevant due to corner cutting
			{
				continue; // skipping this
			}

			double curr_alternative_log_prob = (curr_prev_P_log + curr_Rnm_total_log_prob);

			log_probs.push_back(curr_alternative_log_prob); // add current log prob

			// update max - these will be used in the sum of probs computation:
			if (max_obs_log_prob < curr_alternative_log_prob)
			{
				max_obs_log_prob = curr_alternative_log_prob;
			}
		}
	}

	// compute sum of probs:
	// the idea here is to reduce max_obs_log_prob from all log-probs
	// https://en.wikipedia.org/wiki/LogSumExp
	double threshold_for_exp = std::numeric_limits<double>::min_exponent; // below this there is a risk of an underflow...
	double sum_of_exp_diffs = 0;
	for (size_t i = 0; i < log_probs.size(); i++)
	{
		double curr_log_diff = log_probs[i] - max_obs_log_prob;
		double curr_exp_of_log_diff;
		if (curr_log_diff < threshold_for_exp)
		{
			curr_exp_of_log_diff = exp(threshold_for_exp);
			cout << "when computing the sum of probs we take the exponenet of differences of each log-prob from the max log-prob. In this case the difference was smaller than " << threshold_for_exp << " so we took " << threshold_for_exp << " to avoid an underflow" << endl;
		}
		else
		{
			curr_exp_of_log_diff = exp(curr_log_diff);
		}
		sum_of_exp_diffs += curr_exp_of_log_diff;
	}
	_conditional_alignment_log_total_prob = max_obs_log_prob + log(sum_of_exp_diffs);

	return (_conditional_alignment_log_total_prob);
}

void compute_alignment_dp::fill_S_table_corner_cutting(bool should_take_max, int band_width)
{
	// _S_table contains elements where the first is a chop and the second is a probability.
	// the probability refers to the sampled alignment ending with a match of i and j
	// the chop is the chop that got us there from the previous step
	// it is analogous to _P_table

	vector<size_t> anc_seq = _coded_seqs[0];
	vector<size_t> des_seq = _coded_seqs[1];

	// allocate the table - then we can use the indices
	vector<pair<chop_prob, double>> dummy_row;
	for (size_t col_ind = 0; col_ind < _length_of_des; col_ind++)
	{
		vector<size_t> dummy_des_chars;
		chop_prob dummy_chop(_chop_tables_obj, 'B', 0, 0, 1, 1, dummy_des_chars, _branch_length_t, _is_jc, _quick_jtt_obj);
		dummy_row.push_back(make_pair(dummy_chop, 1.0)); // the value of 1.0 will serve in sanity checks as this is not a legal value for log probability
	}
	for (size_t row_ind = 0; row_ind < _length_of_anc; row_ind++)
	{
		_S_table.push_back(dummy_row);
	}

	// this would allow us to iterate only on chops that exist in the map
	bool should_use_N_chops_accelaration = true;
	map<pair<size_t, size_t>, double> N_chops_map = _chop_tables_obj.get_N_chops_map();
	size_t num_N_elements = N_chops_map.size();
	size_t app_complexity_N_chops = _length_of_anc * _length_of_des * num_N_elements;
	size_t app_complexity_pairs = 0;
	if (band_width < 0)
	{
		app_complexity_pairs = _length_of_anc * _length_of_des * _length_of_anc * _length_of_des;
	}
	else
	{
		size_t width_with_band = abs(int(_length_of_anc) - (int)_length_of_des) + 2 * abs(band_width);
		app_complexity_pairs = (_length_of_anc * width_with_band) * (_length_of_des * width_with_band);
	}
	if (app_complexity_pairs < app_complexity_N_chops)
	{
		should_use_N_chops_accelaration = false;
	}
	cout << "band width parameter is: " << band_width << endl;
	cout << "curr table has " << N_chops_map.size() << " elements in its N table." << endl;
	cout << "estimated number of operations under band-accelartion is: " << app_complexity_pairs << endl;
	cout << "estimated number of operations under chop-accelartion is: " << app_complexity_N_chops << endl;
	cout << "will chop-accelartion be used: " << should_use_N_chops_accelaration << endl;

	vector<vector<size_t>> i_and_j_combs;
	get_combs(band_width, i_and_j_combs);

	for (size_t i_and_j_comb_ind = 0; i_and_j_comb_ind < i_and_j_combs.size(); i_and_j_comb_ind++)
	{
		size_t anc_ind_i = i_and_j_combs[i_and_j_comb_ind][0];
		size_t des_ind_j = i_and_j_combs[i_and_j_comb_ind][1];
		
		size_t anc_matched_char = anc_seq[anc_ind_i]; // match on the right - this is part of the event
		size_t des_matched_char = des_seq[des_ind_j]; // match on the right - this is part of the event

		vector<size_t> inserted_des_chars;
		for (size_t k = 0; k < des_ind_j; k++) // in case des_ind_j = 0, this will be empty
		{
			inserted_des_chars.push_back(des_seq[k]);
		}

		chop_prob curr_Lij_chop(_chop_tables_obj, 'L', anc_ind_i, des_ind_j, anc_matched_char, des_matched_char, inserted_des_chars, _branch_length_t, _is_jc, _quick_jtt_obj);
		double curr_Lij_total_log_prob = curr_Lij_chop.get_chop_log_prob();

		chop_prob curr_path_max_chop = curr_Lij_chop; // can be used for max, if needed
		double curr_iteration_path_max_log_prob = curr_Lij_total_log_prob; // can be used for max, if needed

		vector<chop_prob> chop_options; // can be used for sampling, if needed
		vector<double> log_probs; // can be used for sampling, if needed
		chop_options.push_back(curr_Lij_chop); // can be used for sampling, if needed
		log_probs.push_back(curr_Lij_total_log_prob); // can be used for sampling, if needed

		if ((anc_ind_i > 0) && (des_ind_j > 0))
		{
			// more efficiant to iterate over pairs:
			if (!should_use_N_chops_accelaration)
			{
				vector<vector<size_t>> curr_comb_prev_combs;
				get_prev_combs_for_pair(band_width, curr_comb_prev_combs, anc_ind_i, des_ind_j);

				for (size_t prev_i_j_combs_ind = 0; prev_i_j_combs_ind < curr_comb_prev_combs.size(); prev_i_j_combs_ind++)
				{
					size_t curr_anc_prev_ind = curr_comb_prev_combs[prev_i_j_combs_ind][0];
					size_t curr_des_prev_ind = curr_comb_prev_combs[prev_i_j_combs_ind][1];

					size_t n = anc_ind_i - curr_anc_prev_ind - 1;
					size_t m = des_ind_j - curr_des_prev_ind - 1;

					vector<size_t> inserted_des_chars_internal; // in case m = 0, this will be empty
					if (m > 0)
					{
						// note: (des_ind_j - m) = curr_des_prev_ind + 1
						for (size_t k = (des_ind_j - m); k < des_ind_j; k++)
						{
							inserted_des_chars_internal.push_back(des_seq[k]);
						}
					}
					chop_prob curr_Nnm_chop(_chop_tables_obj, 'N', n, m, anc_matched_char, des_matched_char, inserted_des_chars_internal, _branch_length_t, _is_jc, _quick_jtt_obj);
					double curr_Nnm_total_log_prob = curr_Nnm_chop.get_chop_log_prob();
					double curr_prev_path_S = (_S_table[curr_anc_prev_ind][curr_des_prev_ind]).second; // guaranteed to be computed

					if (curr_prev_path_S > 0) // sanity check - preceding elements should be already computed
					{
						cout << "We have a bug: " << curr_prev_path_S << endl;
					}

					double curr_alternative_log_prob = (curr_prev_path_S + curr_Nnm_total_log_prob);

					chop_options.push_back(curr_Nnm_chop); // can be used for sampling, if needed
					log_probs.push_back(curr_alternative_log_prob); // can be used for sampling, if needed

					if (curr_iteration_path_max_log_prob < curr_alternative_log_prob)
					{
						curr_iteration_path_max_log_prob = curr_alternative_log_prob; // can be used for max, if needed
						curr_path_max_chop = curr_Nnm_chop; // can be used for max, if needed
					}
				}
			}
			else // more efficient to iterate over N_chops_map:
			{
				for (auto const& possible_N_chop : N_chops_map)
				{
					pair<size_t, size_t> curr_ij = possible_N_chop.first;
					size_t n = curr_ij.first;
					size_t m = curr_ij.second;

					int possible_curr_anc_prev_ind = (int)anc_ind_i - (int)n - 1;
					int possible_curr_des_prev_ind = (int)des_ind_j - (int)m - 1;

					if ((possible_curr_anc_prev_ind < 0) || (possible_curr_des_prev_ind < 0))
					{
						continue; // out of bound
					}
					if ((possible_curr_anc_prev_ind >= (int)anc_ind_i) || (possible_curr_des_prev_ind >= (int)des_ind_j))
					{
						continue; // out of bound
					}

					size_t curr_anc_prev_ind = (size_t)possible_curr_anc_prev_ind;
					size_t curr_des_prev_ind = (size_t)possible_curr_des_prev_ind;
					double curr_prev_path_S = (_S_table[curr_anc_prev_ind][curr_des_prev_ind]).second; // computed unless not in band
					if ((curr_prev_path_S > 0) && (band_width == -1)) // sanity check - if no band - all preceding elements should be already computed
					{
						cout << "We have a bug: " << curr_prev_path_S << endl;
					}
					else if (curr_prev_path_S > 0) // not computed, if no bug - this is due to band
					{
						continue; // no need to consider this pair - it's out of band
					}

					// now - the actual computation:
					vector<size_t> inserted_des_chars_internal; // in case m = 0, this will be empty
					if (m > 0)
					{
						// note: (des_ind_j - m) = curr_des_prev_ind + 1
						for (size_t k = (des_ind_j - m); k < des_ind_j; k++)
						{
							inserted_des_chars_internal.push_back(des_seq[k]);
						}
					}
					chop_prob curr_Nnm_chop(_chop_tables_obj, 'N', n, m, anc_matched_char, des_matched_char, inserted_des_chars_internal, _branch_length_t, _is_jc, _quick_jtt_obj);
					double curr_Nnm_total_log_prob = curr_Nnm_chop.get_chop_log_prob();

					double curr_alternative_log_prob = (curr_prev_path_S + curr_Nnm_total_log_prob);

					chop_options.push_back(curr_Nnm_chop); // can be used for sampling, if needed
					log_probs.push_back(curr_alternative_log_prob); // can be used for sampling, if needed

					if (curr_iteration_path_max_log_prob < curr_alternative_log_prob)
					{
						curr_iteration_path_max_log_prob = curr_alternative_log_prob; // can be used for max, if needed
						curr_path_max_chop = curr_Nnm_chop; // can be used for max, if needed
					}
				}
			}
		}

		if (should_take_max) // no sampling, take max
		{
			_S_table[anc_ind_i][des_ind_j] = make_pair(curr_path_max_chop, curr_iteration_path_max_log_prob);
		}
		else // sample
		{
			size_t sampled_ind = sample_index_gumbel_max(log_probs);
			_S_table[anc_ind_i][des_ind_j] = make_pair(chop_options[sampled_ind], log_probs[sampled_ind]);
			//cout << "i:" << anc_ind_i << " j:" << des_ind_j << ", sampled index " << sampled_ind << " out of " << log_probs.size() << " options" << endl;
		}
	}
}

void compute_alignment_dp::fill_P_table_corner_cutting(int band_width)
{
	// _P_table contains the "forward" probabilities (eq 14 of Miklos et al 2004 with zero indices)
	// we work in log-space but sum probabilites. This is possible thanks to:
	// 1. collect into a vector all log - probs(denote each ai)
	// 2. https://en.wikipedia.org/wiki/LogSumExp

	vector<size_t> anc_seq = _coded_seqs[0];
	vector<size_t> des_seq = _coded_seqs[1];

	// allocate the table - then we can use the indices
	_P_table.resize(_length_of_anc);
	for (size_t row_ind = 0; row_ind < _length_of_anc; row_ind++)
	{
		_P_table[row_ind].resize(_length_of_des);
	}
	// initialize with illegitimate values:
	for (size_t i = 0; i < _length_of_anc; i++)
	{
		for (size_t j = 0; j < _length_of_des; j++)
		{
			_P_table[i][j] = 1.0; // not a legitimate log value of probabilities
		}
	}

	// this would allow us to iterate only on chops that exist in the map
	bool should_use_N_chops_accelaration = true;
	map<pair<size_t, size_t>, double> N_chops_map = _chop_tables_obj.get_N_chops_map();
	size_t num_N_elements = N_chops_map.size();
	size_t app_complexity_N_chops = _length_of_anc * _length_of_des * num_N_elements;
	size_t app_complexity_pairs = 0;
	if (band_width < 0)
	{
		app_complexity_pairs = _length_of_anc * _length_of_des * _length_of_anc * _length_of_des;
	}
	else
	{
		app_complexity_pairs = (_length_of_anc * abs(int(_length_of_anc) - (int)_length_of_des) * abs(band_width)) * (_length_of_des * abs(int(_length_of_anc) - (int)_length_of_des) * abs(band_width));
	}
	if (app_complexity_pairs < app_complexity_N_chops)
	{
		should_use_N_chops_accelaration = false;
	}
	cout << "curr table has " << N_chops_map.size() << " elements in its N table." << endl;
	cout << "estimated number of operations under band-accelartion is: " << app_complexity_pairs << endl;
	cout << "estimated number of operations under chop-accelartion is: " << app_complexity_N_chops << endl;
	cout << "will chop-accelartion be used: " << should_use_N_chops_accelaration << endl;

	vector<vector<size_t>> i_and_j_combs;
	get_combs(band_width, i_and_j_combs);

	for (size_t i_and_j_comb_ind = 0; i_and_j_comb_ind < i_and_j_combs.size(); i_and_j_comb_ind++)
	{
		size_t anc_ind_i = i_and_j_combs[i_and_j_comb_ind][0];
		size_t des_ind_j = i_and_j_combs[i_and_j_comb_ind][1];
		
		size_t anc_matched_char = anc_seq[anc_ind_i]; // match on the right - this is part of the event
		size_t des_matched_char = des_seq[des_ind_j]; // match on the right - this is part of the event

		vector<size_t> inserted_des_chars;
		for (size_t k = 0; k < des_ind_j; k++) // in case des_ind_j = 0, this will be empty
		{
			inserted_des_chars.push_back(des_seq[k]);
		}

		chop_prob curr_Lij_chop(_chop_tables_obj, 'L', anc_ind_i, des_ind_j, anc_matched_char, des_matched_char, inserted_des_chars, _branch_length_t, _is_jc, _quick_jtt_obj);
		double curr_Lij_total_log_prob = curr_Lij_chop.get_chop_log_prob();

		vector<double> log_probs; // will be used for summing probs
		log_probs.push_back(curr_Lij_total_log_prob); // add the first log prob
		double max_obs_log_prob = curr_Lij_total_log_prob; // initialize

		if ((anc_ind_i > 0) && (des_ind_j > 0))
		{
			// pairs accleration is more efficiant to iterate over pairs:
			if (!should_use_N_chops_accelaration)
			{
				vector<vector<size_t>> curr_comb_prev_combs;
				get_prev_combs_for_pair(band_width, curr_comb_prev_combs, anc_ind_i, des_ind_j);

				for (size_t prev_i_j_combs_ind = 0; prev_i_j_combs_ind < curr_comb_prev_combs.size(); prev_i_j_combs_ind++)
				{
					size_t curr_anc_prev_ind = curr_comb_prev_combs[prev_i_j_combs_ind][0];
					size_t curr_des_prev_ind = curr_comb_prev_combs[prev_i_j_combs_ind][1];

					size_t n = anc_ind_i - curr_anc_prev_ind - 1;
					size_t m = des_ind_j - curr_des_prev_ind - 1;

					vector<size_t> inserted_des_chars_internal; // in case m = 0, this will be empty
					if (m > 0)
					{
						// note: (des_ind_j - m) = curr_des_prev_ind + 1
						for (size_t k = (des_ind_j - m); k < des_ind_j; k++)
						{
							inserted_des_chars_internal.push_back(des_seq[k]);
						}
					}
					chop_prob curr_Nnm_chop(_chop_tables_obj, 'N', n, m, anc_matched_char, des_matched_char, inserted_des_chars_internal, _branch_length_t, _is_jc, _quick_jtt_obj);
					double curr_Nnm_total_log_prob = curr_Nnm_chop.get_chop_log_prob();
					double curr_prev_P_log = (_P_table[curr_anc_prev_ind][curr_des_prev_ind]); // guaranteed to be computed

					if (curr_prev_P_log > 0) // sanity check - preceding elements should be already computed
					{
						cout << "We have a bug: " << curr_prev_P_log << endl;
					}

					double curr_alternative_log_prob = (curr_prev_P_log + curr_Nnm_total_log_prob);

					log_probs.push_back(curr_alternative_log_prob); // add current log prob

					// update max - these will be used in the sum of probs computation:
					if (max_obs_log_prob < curr_alternative_log_prob)
					{
						max_obs_log_prob = curr_alternative_log_prob;
					}
				}
			}
			else // more efficiant to iterate over N_chops_map:
			{
				for (auto const& possible_N_chop : N_chops_map)
				{
					pair<size_t, size_t> curr_ij = possible_N_chop.first;
					size_t n = curr_ij.first;
					size_t m = curr_ij.second;

					int possible_curr_anc_prev_ind = (int)anc_ind_i - (int)n - 1;
					int possible_curr_des_prev_ind = (int)des_ind_j - (int)m - 1;

					if ((possible_curr_anc_prev_ind < 0) || (possible_curr_des_prev_ind < 0))
					{
						continue; // out of bound
					}
					if ((possible_curr_anc_prev_ind >= (int)anc_ind_i) || (possible_curr_des_prev_ind >= (int)des_ind_j))
					{
						continue; // out of bound
					}

					size_t curr_anc_prev_ind = (size_t)possible_curr_anc_prev_ind;
					size_t curr_des_prev_ind = (size_t)possible_curr_des_prev_ind;
					double curr_prev_P_log = (_P_table[curr_anc_prev_ind][curr_des_prev_ind]); // computed unless not in band
					if ((curr_prev_P_log > 0) && (band_width == -1)) // sanity check - if no band - all preceding elements should be already computed
					{
						cerr << "We have a bug: " << curr_prev_P_log << endl;
					}
					else if (curr_prev_P_log > 0) // not computed, if no bug - this is due to band
					{
						continue; // no need to consider this pair - it's out of band
					}

					// now - the actual computation:
					vector<size_t> inserted_des_chars_internal; // in case m = 0, this will be empty
					if (m > 0)
					{
						// note: (des_ind_j - m) = curr_des_prev_ind + 1
						for (size_t k = (des_ind_j - m); k < des_ind_j; k++)
						{
							inserted_des_chars_internal.push_back(des_seq[k]);
						}
					}
					chop_prob curr_Nnm_chop(_chop_tables_obj, 'N', n, m, anc_matched_char, des_matched_char, inserted_des_chars_internal, _branch_length_t, _is_jc, _quick_jtt_obj);
					double curr_Nnm_total_log_prob = curr_Nnm_chop.get_chop_log_prob();

					double curr_alternative_log_prob = (curr_prev_P_log + curr_Nnm_total_log_prob);

					log_probs.push_back(curr_alternative_log_prob); // add current log prob

					// update max - these will be used in the sum of probs computation:
					if (max_obs_log_prob < curr_alternative_log_prob)
					{
						max_obs_log_prob = curr_alternative_log_prob;
					}
				}
			}
		}

		// compute sum of probs:
		// the idea here is to reduce max_obs_log_prob from all log-probs
		// we can then do exp(a - max_obs_log_prob), where 'a' is each of the log-probs
		// then we take log(exp(a0 - max_obs_log_prob) + exp(a1 - max_obs_log_prob) + ... )
		// https://en.wikipedia.org/wiki/LogSumExp
		double threshold_for_exp = std::numeric_limits<double>::min_exponent; // below this there is a risk of an underflow...
		double sum_of_exp_reduced_log_probs = 0.0;
		for (size_t i = 0; i < log_probs.size(); i++)
		{
			//cout << "curr log-prob is: " << log_probs[i] << " and max log-prob is: " << max_obs_log_prob << endl;
			double curr_log_diff = log_probs[i] - max_obs_log_prob;
			double curr_exp_of_log_diff;
			if (curr_log_diff < threshold_for_exp)
			{
				curr_exp_of_log_diff = exp(threshold_for_exp);
				cout << "curr log-prob is: " << log_probs[i] << " and max log-prob is: " << max_obs_log_prob << endl;
				cout << "curr anc_ind_i is: " << anc_ind_i << " and des_ind_j is: " << des_ind_j << endl;
				cout << "when computing the sum of probs we take the exponenet of differences of each log-prob from the max log-prob. In this case the difference was smaller than " << threshold_for_exp << " so we took " << threshold_for_exp << " to avoid an underflow" << endl;
			}
			else
			{
				curr_exp_of_log_diff = exp(curr_log_diff);
			}
			sum_of_exp_reduced_log_probs += curr_exp_of_log_diff;
		}
		_P_table[anc_ind_i][des_ind_j] = log(sum_of_exp_reduced_log_probs) + max_obs_log_prob;
	}

	_is_P_computed = true; // now it is computed - no need to recompute in the future
}

void compute_alignment_dp::fill_X_table_corner_cutting(int band_width)
{
	// _X_table contains the "backward" probabilities (analog of eq 14 of Miklos et al 2004 with zero indices)
	// we work in log-space but sum probabilites. This is possible thanks to:
	// 1. collect into a vector all log - probs(denote each ai)
	// https://en.wikipedia.org/wiki/LogSumExp

	vector<size_t> anc_seq = _coded_seqs[0];
	vector<size_t> des_seq = _coded_seqs[1];

	// allocate the table - then we can use the indices
	_X_table.resize(_length_of_anc);
	for (size_t row_ind = 0; row_ind < _length_of_anc; row_ind++)
	{
		_X_table[row_ind].resize(_length_of_des);
	}
	// initialize with illegitimate values:
	for (size_t i = 0; i < _length_of_anc; i++)
	{
		for (size_t j = 0; j < _length_of_des; j++)
		{
			_X_table[i][j] = 1.0; // not a legitimate log value of probabilities
		}
	}

	// this would allow us to iterate only on chops that exist in the map
	bool should_use_N_chops_accelaration = true;
	map<pair<size_t, size_t>, double> N_chops_map = _chop_tables_obj.get_N_chops_map();
	size_t num_N_elements = N_chops_map.size();
	size_t app_complexity_N_chops = _length_of_anc * _length_of_des * num_N_elements;
	size_t app_complexity_pairs = 0;
	if (band_width < 0)
	{
		app_complexity_pairs = _length_of_anc * _length_of_des * _length_of_anc * _length_of_des;
	}
	else
	{
		app_complexity_pairs = (_length_of_anc * abs(int(_length_of_anc) - (int)_length_of_des) * abs(band_width)) * (_length_of_des * abs(int(_length_of_anc) - (int)_length_of_des) * abs(band_width));
	}
	if (app_complexity_pairs < app_complexity_N_chops)
	{
		should_use_N_chops_accelaration = false;
	}
	cout << "curr table has " << N_chops_map.size() << " elements in its N table." << endl;
	cout << "estimated number of operations under band-accelartion is: " << app_complexity_pairs << endl;
	cout << "estimated number of operations under chop-accelartion is: " << app_complexity_N_chops << endl;
	cout << "will chop-accelartion be used: " << should_use_N_chops_accelaration << endl;

	vector<vector<size_t>> i_and_j_combs_backward;
	get_backward_combs(band_width, i_and_j_combs_backward);

	for (size_t i_and_j_comb_ind = 0; i_and_j_comb_ind < i_and_j_combs_backward.size(); i_and_j_comb_ind++)
	{
		size_t anc_ind_i = i_and_j_combs_backward[i_and_j_comb_ind][0];
		size_t des_ind_j = i_and_j_combs_backward[i_and_j_comb_ind][1];

		vector<size_t> inserted_des_chars; // in case (_length_of_des - des_ind_j - 1) = 0, this will be empty
		if ((_length_of_des - des_ind_j - 1) > 0)
		{
			for (size_t k = (des_ind_j + 1); k < _length_of_des; k++)
			{
				inserted_des_chars.push_back(des_seq[k]);
			}
		}

		chop_prob curr_R_chop(_chop_tables_obj, 'R', (_length_of_anc - anc_ind_i - 1), (_length_of_des - des_ind_j - 1), 0, 0, inserted_des_chars, _branch_length_t, _is_jc, _quick_jtt_obj);
		double curr_R_total_log_prob = curr_R_chop.get_chop_log_prob();

		vector<double> log_probs; // will be used for summing probs
		log_probs.push_back(curr_R_total_log_prob); // add the first log prob
		double max_obs_log_prob = curr_R_total_log_prob; // initialize

		if ((anc_ind_i < (_length_of_anc - 1)) && (des_ind_j < (_length_of_des - 1)))
		{
			// pairs accleration is more efficiant to iterate over pairs:
			if (!should_use_N_chops_accelaration)
			{
				vector<vector<size_t>> curr_comb_prev_combs_backward;
				get_backward_prev_combs_for_pair(band_width, curr_comb_prev_combs_backward, anc_ind_i, des_ind_j);

				for (size_t prev_i_j_combs_ind = 0; prev_i_j_combs_ind < curr_comb_prev_combs_backward.size(); prev_i_j_combs_ind++)
				{
					size_t curr_anc_prev_ind = curr_comb_prev_combs_backward[prev_i_j_combs_ind][0];
					size_t curr_des_prev_ind = curr_comb_prev_combs_backward[prev_i_j_combs_ind][1];

					size_t num_del_from_anc = curr_anc_prev_ind - anc_ind_i + 1; // always >= 0
					size_t num_ins_in_des = curr_des_prev_ind - des_ind_j + 1; // always >= 0

					vector<size_t> inserted_des_chars_internal; // in case num_ins_in_des = 0, this will be empty
					if (num_ins_in_des > 0)
					{
						for (size_t k = (des_ind_j + 1); k < num_ins_in_des; k++)
						{
							inserted_des_chars_internal.push_back(des_seq[k]);
						}
					}

					size_t anc_matched_char = anc_seq[curr_anc_prev_ind];
					size_t des_matched_char = des_seq[curr_des_prev_ind];

					chop_prob curr_Nnm_chop(_chop_tables_obj, 'N', num_del_from_anc, num_ins_in_des, anc_matched_char, des_matched_char, inserted_des_chars_internal, _branch_length_t, _is_jc, _quick_jtt_obj);
					double curr_Nnm_total_log_prob = curr_Nnm_chop.get_chop_log_prob();
					double curr_prev_X_log = _X_table[curr_anc_prev_ind][curr_des_prev_ind]; // guaranteed to be computed

					if (curr_prev_X_log > 0) // sanity check - preceding elements should be already computed
					{
						cout << "We have a bug: " << curr_prev_X_log << endl;
					}

					double curr_alternative_log_prob = (curr_prev_X_log + curr_Nnm_total_log_prob);

					log_probs.push_back(curr_alternative_log_prob); // add current log prob

																	// update max - these will be used in the sum of probs computation:
					if (max_obs_log_prob < curr_alternative_log_prob)
					{
						max_obs_log_prob = curr_alternative_log_prob;
					}
				}
			}
			else // more efficient to iterate over N_chops_map
			{
				for (auto const& possible_N_chop : N_chops_map)
				{
					pair<size_t, size_t> curr_ij = possible_N_chop.first;
					size_t num_del_from_anc = curr_ij.first;
					size_t num_ins_in_des = curr_ij.second;

					int possible_curr_anc_prev_ind = (int)anc_ind_i + (int)num_del_from_anc - 1;
					int possible_curr_des_prev_ind = (int)des_ind_j + (int)num_ins_in_des - 1;

					if ((possible_curr_anc_prev_ind > ((int)_length_of_anc - 1)) || (possible_curr_des_prev_ind > ((int)_length_of_des - 1)))
					{
						continue; // out of bound
					}
					if ((possible_curr_anc_prev_ind < ((int)anc_ind_i + 1)) || (possible_curr_des_prev_ind < ((int)des_ind_j + 1)))
					{
						continue; // out of bound
					}

					size_t curr_anc_prev_ind = (size_t)possible_curr_anc_prev_ind;
					size_t curr_des_prev_ind = (size_t)possible_curr_des_prev_ind;
					double curr_prev_X_log = _X_table[curr_anc_prev_ind][curr_des_prev_ind]; // computed unless not in band
					if ((curr_prev_X_log > 0) && (band_width == -1)) // sanity check - if no band - all preceding elements should be already computed
					{
						cout << "We have a bug: " << curr_prev_X_log << endl;
					}
					else if (curr_prev_X_log > 0) // not computed, if no bug - this is due to band
					{
						continue; // no need to consider this pair - it's out of band
					}

					// now - the actual computation:
					vector<size_t> inserted_des_chars_internal; // in case num_ins_in_des = 0, this will be empty
					if (num_ins_in_des > 0)
					{
						for (size_t k = (des_ind_j + 1); k < num_ins_in_des; k++)
						{
							inserted_des_chars_internal.push_back(des_seq[k]);
						}
					}

					size_t anc_matched_char = anc_seq[curr_anc_prev_ind];
					size_t des_matched_char = des_seq[curr_des_prev_ind];

					chop_prob curr_Nnm_chop(_chop_tables_obj, 'N', num_del_from_anc, num_ins_in_des, anc_matched_char, des_matched_char, inserted_des_chars_internal, _branch_length_t, _is_jc, _quick_jtt_obj);
					double curr_Nnm_total_log_prob = curr_Nnm_chop.get_chop_log_prob();

					double curr_alternative_log_prob = (curr_prev_X_log + curr_Nnm_total_log_prob);

					log_probs.push_back(curr_alternative_log_prob); // add current log prob

					// update max - these will be used in the sum of probs computation:
					if (max_obs_log_prob < curr_alternative_log_prob)
					{
						max_obs_log_prob = curr_alternative_log_prob;
					}
				}
			}		
		}

		// compute sum of probs:
		// the idea here is to reduce max_obs_log_prob from all log-probs
		// we can then do exp(a - max_obs_log_prob), where 'a' is each of the log-probs
		// then we take log(exp(a0 - max_obs_log_prob) + exp(a1 - max_obs_log_prob) + ... )
		// https://en.wikipedia.org/wiki/LogSumExp
		double threshold_for_exp = std::numeric_limits<double>::min_exponent; // below this there is a risk of an underflow...
		double sum_of_exp_reduced_log_probs = 0.0;
		for (size_t i = 0; i < log_probs.size(); i++)
		{
			double curr_log_diff = log_probs[i] - max_obs_log_prob;
			double curr_exp_of_log_diff;
			if (curr_log_diff < threshold_for_exp)
			{
				curr_exp_of_log_diff = exp(threshold_for_exp);
				cout << "when computing the sum of probs we take the exponenet of differences of each log-prob from the max log-prob. In this case the difference was smaller than " << threshold_for_exp << " so we took " << threshold_for_exp << " to avoid an underflow" << endl;
			}
			else
			{
				curr_exp_of_log_diff = exp(curr_log_diff);
			}
			sum_of_exp_reduced_log_probs += curr_exp_of_log_diff;
		}
		_X_table[anc_ind_i][des_ind_j] = log(sum_of_exp_reduced_log_probs) + max_obs_log_prob;
	}

	_is_X_computed = true; // now it is computed - no need to recompute in the future
}