// code by Eli Levy Karin

#include "chop_prob.h"

chop_prob & chop_prob::operator=(const chop_prob & other)
{
	_quick_jtt_obj = other._quick_jtt_obj;
	_chop_type = other._chop_type;
	_branch_length_t = other._branch_length_t;
	_i = other._i;
	_j = other._j;
	_des_inserted_chars = other._des_inserted_chars;
	_anc_matched_char = other._anc_matched_char;
	_des_matched_char = other._des_matched_char;
	_chop_log_prob = other._chop_log_prob;
	_sub_log_prob = other._sub_log_prob;
	_indel_log_prob = other._indel_log_prob;
	return *this;
}

void chop_prob::print_chop_info() const
{
	cout << "current chop has the following features:" << endl;
	cout << "type = " << _chop_type << endl;
	cout << "i = " << _i << endl;
	cout << "j = " << _j << endl;
	if (_j != 0)
	{
		cout << "inserted des chars are: ";
		for (size_t anc_ins_char_ind = 0; anc_ins_char_ind < _des_inserted_chars.size(); anc_ins_char_ind++)
		{
			cout << _des_inserted_chars[anc_ins_char_ind];
			if ((anc_ins_char_ind + 1) != _des_inserted_chars.size())
			{
				cout << ",";
			}
		}
		cout << endl;
	}
	cout << "matched char anc = " << _anc_matched_char << endl;
	cout << "matched char des = " << _des_matched_char << endl;
	cout << "substitution log prob = " << _sub_log_prob << endl;
	cout << "indel log prob = " << _indel_log_prob << endl;
	cout << "total log prob = " << _chop_log_prob << endl;
	cout << endl;
}

void chop_prob::compute_JC_substitution_prob()
{
	size_t num_inserted_chars_des = _des_inserted_chars.size();
	if (_j != num_inserted_chars_des)
	{
		cerr << "j: " << _j << " and num_inserted_chars_des: " << num_inserted_chars_des << " should be equal!" << endl;
		exit(1);
	}
	if (_j > 0) // there was at least one inserted character in the descendant
	{
		for (size_t ins_char_ind = 0; ins_char_ind < _des_inserted_chars.size(); ins_char_ind++)
		{
			_sub_log_prob += log(0.25);
		}
	}
	if ((_chop_type == 'N') || (_chop_type == 'L')) // only L and N chops take a match into account
	{
		if ((_anc_matched_char == 0) || (_des_matched_char == 0)) // this should not happen - sanity check!
		{
			cerr << "something is off! a matched char is 0, which denotes a gap while the chop type is " << _chop_type << endl;
			exit(1);
		}
		double trans_prob = 0.25;
		if (_anc_matched_char == _des_matched_char)
		{
			trans_prob += 0.75 * exp((-4 * _branch_length_t / 3));
		}
		else
		{
			trans_prob -= 0.25 * exp((-4 * _branch_length_t / 3));
		}
		_sub_log_prob += log(trans_prob);
	}
}

void chop_prob::compute_JTT_substitution_prob()
{
	size_t num_inserted_chars_des = _des_inserted_chars.size();
	if (_j != num_inserted_chars_des)
	{
		cerr << "j: " << _j << " and num_inserted_chars_des: " << num_inserted_chars_des << " should be equal!" << endl;
		exit(1);
	}
	if (_j > 0) // there was at least one inserted character in the descendant
	{
		for (size_t ins_char_ind = 0; ins_char_ind < _des_inserted_chars.size(); ins_char_ind++)
		{
			double aa_freq = _quick_jtt_obj.freq(_des_inserted_chars[ins_char_ind]);
			_sub_log_prob += log(aa_freq);
		}
	}
	if ((_chop_type == 'N') || (_chop_type == 'L')) // only L and N chops take a match into account
	{
		if ((_anc_matched_char == 0) || (_des_matched_char == 0)) // this should not happen - sanity check!
		{
			cerr << "something is off! a matched char is 0, which denotes a gap while the chop type is " << _chop_type << endl;
			exit(1);
		}
		double trans_prob = _quick_jtt_obj.get_pij_t(_anc_matched_char, _des_matched_char, _branch_length_t);
		_sub_log_prob += log(trans_prob);
	}
}

void chop_prob::compute_chop_prob()
{
	if (_is_jc)
	{
		compute_JC_substitution_prob(); // will compute _sub_prob, DNA
	}
	else
	{
		compute_JTT_substitution_prob(); // will compute _sub_prob, AA
	}
	
	_indel_log_prob = log(_chop_tables_obj.get_chop_prob(_chop_type, _i, _j));
	_chop_log_prob = _sub_log_prob + _indel_log_prob;
}
