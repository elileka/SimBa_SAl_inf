// code by Eli Levy Karin

#include "alignment.h"

double alignment::get_alignment_log_probability_cond_on_anc() const
{
	double alignment_log_prob = 0;
	for (size_t chop_ind = 0; chop_ind < _aligment_chops.size(); chop_ind++)
	{
		//_aligment_chops[chop_ind].print_chop_info(); // for debug
		double curr_chop_prob = _aligment_chops[chop_ind].get_chop_log_prob();
		alignment_log_prob += curr_chop_prob;
	}
	return alignment_log_prob;
}

void alignment::break_alignment_to_chops()
{
	bool leftmost = true;
	bool saw_a_match = false;
	size_t curr_i = 0;
	size_t curr_j = 0;
	vector<size_t> inserted_des_chars;
	size_t anc_matched_char = 0; // if chop is of type L or N, we change this to the actual match
	size_t des_matched_char = 0; // if chop is of type L or N, we change this to the actual match

	for (size_t pos_ind = 0; pos_ind < _coded_alignment.size(); pos_ind++)
	{
		if ((_coded_alignment[pos_ind][0] != 0) && (_coded_alignment[pos_ind][1] != 0)) // match
		{
			saw_a_match = true;
			anc_matched_char = _coded_alignment[pos_ind][0];
			des_matched_char = _coded_alignment[pos_ind][1];
			if (leftmost)
			{
				// L chop
				chop_prob curr_chop(_chop_tables_obj, 'L', curr_i, curr_j, anc_matched_char, des_matched_char, inserted_des_chars, _branch_length_t, _is_jc, _quick_jtt_obj);
				_aligment_chops.push_back(curr_chop);
				leftmost = false;
			}
			else
			{
				// N chop
				chop_prob curr_chop(_chop_tables_obj, 'N', curr_i, curr_j, anc_matched_char, des_matched_char, inserted_des_chars, _branch_length_t, _is_jc, _quick_jtt_obj);
				_aligment_chops.push_back(curr_chop);
			}
			curr_i = 0; // initialize
			curr_j = 0; // initialize
			inserted_des_chars.clear(); // initialize
		}
		else if (_coded_alignment[pos_ind][1] == 0) // deleted from anc
		{
			curr_i++;
		}
		else if (_coded_alignment[pos_ind][0] == 0) // inserted in des
		{
			curr_j++;
			inserted_des_chars.push_back(_coded_alignment[pos_ind][1]);
		}
	}

	if (saw_a_match) // R chop
	{
		chop_prob curr_chop(_chop_tables_obj, 'R', curr_i, curr_j, 0, 0, inserted_des_chars, _branch_length_t, _is_jc, _quick_jtt_obj);
		_aligment_chops.push_back(curr_chop);
	}
	else // B chop
	{
		chop_prob curr_chop(_chop_tables_obj, 'B', curr_i, curr_j, 0, 0, inserted_des_chars, _branch_length_t, _is_jc, _quick_jtt_obj);
		_aligment_chops.push_back(curr_chop);
	}

}
