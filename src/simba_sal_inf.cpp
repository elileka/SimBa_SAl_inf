// code by Eli Levy Karin

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <random>
#include <map>
#include <cmath>
#include <ctime>
#include "read_chop_tables.h"
#include "alignment.h"
#include "compute_alignment_dp.h"
#include "quick_jtt.h"
#include "read_input_seqs.h"

using namespace std;

void print_seqs_fasta(vector<string> & aligned_seqs, ofstream & myfilestream)
{
	string anc_aligned = aligned_seqs[0];
	string des_aligned = aligned_seqs[1];

	myfilestream << ">anc" << endl;
	myfilestream << anc_aligned << endl;
	myfilestream << ">des" << endl;
	myfilestream << des_aligned << endl << endl;
}

vector<double> get_ts_to_consider(double t_MP)
{
	size_t num_ts_to_consider = 10;
	double diff_between_ts_in_table_world = 0.01;

	t_MP = floor(t_MP * 100.00 + 0.5) / 100.00;
	double t_min = t_MP;
	double t_max = 2 * t_min;
	double l_t_min = log(t_min);
	double l_t_max = log(t_max);
	double l_range = l_t_max - l_t_min;
	double l_interval = l_range / (num_ts_to_consider - 1);
	double l_t_next = l_t_min;
	vector<double> l_ts_to_consider;
	while (l_t_next <= l_t_max)
	{
		l_ts_to_consider.push_back(l_t_next);
		l_t_next = l_t_next + l_interval;
	}
	vector<double> ts_to_consider;
	ts_to_consider.push_back(t_min);
	for (size_t i = 1; i < l_ts_to_consider.size(); i++)
	{
		double curr_t = exp(l_ts_to_consider[i]);
		curr_t = floor(curr_t * 100.00 + 0.5) / 100.00;
		if (curr_t - ts_to_consider.back() >= diff_between_ts_in_table_world)
		{
			ts_to_consider.push_back(curr_t);
		}
	}
	return (ts_to_consider);
}

string getCmdOption(int num_args, const char* argv[], const std::string& option)
{
	string val;
	for (int i = 0; i < num_args; ++i)
	{
		string arg = argv[i];
		size_t start_position_of_param_name = arg.find(option); // found is the start position of the match
		if (start_position_of_param_name == 0)
		{
			size_t start_pos_of_value = start_position_of_param_name + option.size();
			val = arg.substr(start_pos_of_value);
			return val;
		}
	}
	return val;
}

int main(int argc, const char * argv[])
{
	if (argc < 2)
	{
		cerr << "The program takes parameter values using '=' pairs. The order does not matter." << endl;
		cerr << "The parameter names are: " << endl;
		cerr << "input_fasta_file (full path to the input sequences can be aligned or unaligned)" << endl;
		cerr << "file_with_paths_to_chop_tables (full path to a file listing files with chop tables)" << endl;
		cerr << "type_seq: DNA/AA (defaults to AA)" << endl;
		cerr << "path_to_chebi_JTT_coef_file (full path, only if AA)" << endl;
		cerr << "type_file: aligned/unaligned (defaults to unaligned)" << endl;
		cerr << "number_to_sample (number of alignemnts to sample, a non-negative int, defaults to 0)" << endl;
		cerr << "band_width (an int for corner cutting, defaults to -1 --> no corner cutting)" << endl;
		cerr << "opt_method: gotoh/full (defaults to gotoh)" << endl;
		cerr << "out_files_prefix" << endl;

		cerr << "Example usage: " << endl << argv[0] << " input_fasta_file=INPUT.fas file_with_paths_to_chop_tables=LIST_OF_TABLE_FILES.txt type_seq=AA path_to_chebi_JTT_coef_file=PATH/chebi_coef.txt type_file=aligned out_files_prefix=OUT_INF_" << endl;
		return 1;
	}

	// collection of user parameter values (some may be missing)
	string input_fasta_file = getCmdOption(argc, argv, "input_fasta_file=");
	string file_with_paths_to_chop_tables = getCmdOption(argc, argv, "file_with_paths_to_chop_tables=");
	string type_seq = getCmdOption(argc, argv, "type_seq=");
	string path_to_chebi_JTT_coef_file = getCmdOption(argc, argv, "path_to_chebi_JTT_coef_file=");
	string type_file = getCmdOption(argc, argv, "type_file=");
	string number_to_sample_str = getCmdOption(argc, argv, "number_to_sample=");
	string band_width_str = getCmdOption(argc, argv, "band_width=");
	string opt_method = getCmdOption(argc, argv, "opt_method=");
	string out_files_prefix = getCmdOption(argc, argv, "out_files_prefix=");

	// special treatment of parameters with default values:
	bool is_jc = false; // defaults to AA
	if (type_seq == "DNA")
	{
		is_jc = true;
	}
	else
	{
		if (file_with_paths_to_chop_tables == "")
		{
			cerr << "AA type seq is chosen but path file_with_paths_to_chop_tables is not provided!" << endl;
			return 1;
		}
	}
	quick_jtt quick_jtt_obj(path_to_chebi_JTT_coef_file);
	if (!is_jc) // sanity checks in case AA with JTT was chosen
	{
		double stationary_prob_of_A = quick_jtt_obj.freq(1); // sanity
		double prob_transition_A_to_R_in_t_05 = quick_jtt_obj.get_pij_t(1, 2, 0.5); // sanity
		cout << "sanity check - JTT object: PI(A) = " << stationary_prob_of_A << endl; // sanity
		cout << "sanity check - JTT object: P(A->R|t = 0.5) = " << prob_transition_A_to_R_in_t_05 << endl; // sanity

		if (!(stationary_prob_of_A > 0))
		{
			cerr << "sanity check failed - JTT object: PI(A) = " << stationary_prob_of_A << endl; // sanity
			exit(1);
		}
		if (!(prob_transition_A_to_R_in_t_05 > 0))
		{
			cerr << "sanity check - JTT object: P(A->R|t = 0.5) = " << prob_transition_A_to_R_in_t_05 << endl; // sanity
			exit(1);
		}
	}
	int band_width = -1; // default band_width = -1 --> no corner cutting
	if (band_width_str != "")
	{
		band_width = atoi(band_width_str.c_str());
	}
	bool should_gotoh_opt = true; // defualts to gotoh-based optimization
	if (opt_method == "full")
	{
		should_gotoh_opt = false;
	}
	bool is_aligned = false; // defualts to unaligned input
	if (type_file == "aligned")
	{
		is_aligned = true;
	}
	int number_alignments_to_sample = 0; // defualts to 0
	if (number_to_sample_str != "")
	{
		number_alignments_to_sample = atoi(number_to_sample_str.c_str());
	}
	// prepare file names prefix
	stringstream ss_pref;
	ss_pref << out_files_prefix;

	// prepare output log file name
	stringstream log_file;
	log_file << ss_pref.str() << "inference_log.txt";
	std::string log_file_string = log_file.str();
	ofstream myLog;
	myLog.open(log_file_string);

	// prepare output result file name
	stringstream result_file;
	result_file << ss_pref.str() << "inference_result.txt";
	std::string result_file_string = result_file.str();

	// initialize
	double conditional_best_alignment_log_prob = 1.0;
	string best_table_best_alignment = "";

	vector<string> best_alignment; // according to the table most fitting for the best alignment
	double best_branch_length_t; // according to the table most fitting for the best alignment

	read_input_seqs input_seqs_obj(input_fasta_file, is_aligned, is_jc);
	vector<string> input_seqs;
	input_seqs.push_back(input_seqs_obj.get_str_orig_seq(0));
	input_seqs.push_back(input_seqs_obj.get_str_orig_seq(1));
	
	// get Gotoh alignment - we will choose the parameters based on it:
	double conditional_gotoh_alignment_log_prob = 1.0;
	double gotoh_branch_length_t;

	string best_table_gotoh_alignment = "";
	vector<vector<size_t>> gotoh_coded_alignment;
	vector<string> gotoh_aligned_seqs;
	vector<double> ts_to_consider;
	if (should_gotoh_opt)
	{
		input_seqs_obj.compute_gotoh_alignment();
		gotoh_coded_alignment = input_seqs_obj.get_gotoh_coded_alignment();
		gotoh_aligned_seqs.push_back(input_seqs_obj.get_str_from_coded_alignment(0, gotoh_coded_alignment));
		gotoh_aligned_seqs.push_back(input_seqs_obj.get_str_from_coded_alignment(1, gotoh_coded_alignment));
		double estimate_branch_length_t_MP = input_seqs_obj.get_estimated_branch_length_MP();
		cout << "gotoh MP estimate of branch length: " << estimate_branch_length_t_MP << endl;
		ts_to_consider = get_ts_to_consider(estimate_branch_length_t_MP);
		cout << "these t values will be considered: " << endl;
		for (size_t i = 0; i < ts_to_consider.size(); i++)
		{
			cout << ts_to_consider[i] << ", ";
		}
		cout << endl;
	}

	// for time measurement:
	size_t number_of_computed_dp_alignments = 0;
	double elapsed_secs_in_all_dp_procedures = 0;
	
	// go over all tables in the tables file - parameter estimation
	ifstream tables_file(file_with_paths_to_chop_tables);
	string line;
	clock_t begin_param_opt = clock(); // start measureing time
	while (getline(tables_file, line))
	{
		stringstream linestream(line);

		string curr_chop_tables_file_string;
		double r_param;
		double basic_mu;
		double basic_gamma;
		size_t max_indel_length;
		double branch_length_t;

		getline(linestream, curr_chop_tables_file_string, '\t'); // read up-to the first tab (discard tab)
		linestream >> r_param >> basic_mu >> basic_gamma >> max_indel_length >> branch_length_t;
		cout << "sanity check - curr table has branch_length_t: " << branch_length_t << endl; // table validation

		read_chop_tables chop_tables_obj(curr_chop_tables_file_string);
		double prob_N_0_0 = chop_tables_obj.get_chop_prob('N', 0, 0); // table validation
		cout << "sanity check - curr table has prob_N_0_0: " << prob_N_0_0 << endl; // table validation

		if (should_gotoh_opt) // optimize with respect to a fixed Gotoh alignment
		{
			double epsilon = 0.0001;
			for (size_t ind = 0; ind < ts_to_consider.size(); ind++)
			{
				if (abs(ts_to_consider[ind] - branch_length_t) < epsilon)
				{
					myLog << "working on table: " << curr_chop_tables_file_string << endl;
					alignment gotoh_alignment(chop_tables_obj, gotoh_coded_alignment, branch_length_t, is_jc, quick_jtt_obj);
					double curr_table_conditional_gotoh_alignment_log_prob = gotoh_alignment.get_alignment_log_probability_cond_on_anc();
					if ((conditional_gotoh_alignment_log_prob > 0) || (curr_table_conditional_gotoh_alignment_log_prob > conditional_gotoh_alignment_log_prob))
					{
						conditional_gotoh_alignment_log_prob = curr_table_conditional_gotoh_alignment_log_prob;
						best_table_gotoh_alignment = curr_chop_tables_file_string;
						gotoh_branch_length_t = branch_length_t;
					}
					myLog << "gotoh alignment LL with current table: " << curr_table_conditional_gotoh_alignment_log_prob << endl;
					myLog << "gotoh alignment fasta: " << endl;
					print_seqs_fasta(gotoh_aligned_seqs, myLog);
				}
			}
		}
		else // optimize by considering all tables in file (each - a dp procedure)
		{
			myLog << "working on table: " << curr_chop_tables_file_string << endl;
			// prepare dp object:
			compute_alignment_dp dp_alignment(chop_tables_obj, input_seqs_obj, branch_length_t, is_jc, quick_jtt_obj);

			// compute best alignment on dp object:
			double curr_best_alignment_log_prob;
			bool should_take_max = true;
			
			clock_t begin = clock(); // start measureing time
			vector<string> curr_best_alignment = dp_alignment.get_sampled_alignment(curr_best_alignment_log_prob, should_take_max, band_width);
			clock_t end = clock(); // end measuring time
			double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
			elapsed_secs_in_all_dp_procedures = elapsed_secs_in_all_dp_procedures + elapsed_secs;
			number_of_computed_dp_alignments++;

			if ((conditional_best_alignment_log_prob > 0) || (curr_best_alignment_log_prob > conditional_best_alignment_log_prob))
			{
				conditional_best_alignment_log_prob = curr_best_alignment_log_prob;
				best_table_best_alignment = curr_chop_tables_file_string;
				best_alignment = curr_best_alignment;
				best_branch_length_t = branch_length_t;
			}
			myLog << "best alignment LL with current table: " << curr_best_alignment_log_prob << endl;
			myLog << "best alignment fasta: " << endl;
			print_seqs_fasta(curr_best_alignment, myLog);

			// 30.01.2018 just for posterior probability computation:
			//double conditional_alignment_total_log_prob = dp_alignment.compute_conditional_alignment_total_log_prob(band_width);
			//myLog << "conditional alignment total LL with current table: " << conditional_alignment_total_log_prob << endl << endl;

		}
	}
	clock_t end_param_opt = clock(); // end measuring time
	double elapsed_secs_param_opt = double(end_param_opt - begin_param_opt) / CLOCKS_PER_SEC;

	// end parameter optimization stage!!!

	if (should_gotoh_opt) // now compute ML-PWA according to gotoh table
	{
		if (conditional_gotoh_alignment_log_prob > 0)
		{
			cerr << "looks like no tables were used with Gotoh - probably t didn't match... check log" << endl;
			exit(1);
		}
		// create chop table object:
		read_chop_tables gotoh_chop_tables_obj(best_table_gotoh_alignment);
		// prepare dp object:
		compute_alignment_dp dp_alignment(gotoh_chop_tables_obj, input_seqs_obj, gotoh_branch_length_t, is_jc, quick_jtt_obj);
		// compute best alignment on dp object:
		bool should_take_max = true;

		clock_t begin = clock(); // start measureing time
		best_alignment = dp_alignment.get_sampled_alignment(conditional_best_alignment_log_prob, should_take_max, band_width);
		clock_t end = clock(); // end measuring time
		double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
		elapsed_secs_in_all_dp_procedures = elapsed_secs_in_all_dp_procedures + elapsed_secs;
		number_of_computed_dp_alignments++;

		best_table_best_alignment = best_table_gotoh_alignment;
		best_branch_length_t = gotoh_branch_length_t;

		myLog << "best alignment LL with gotoh table: " << conditional_best_alignment_log_prob << endl;
		myLog << "best alignment fasta: " << endl;
		print_seqs_fasta(best_alignment, myLog);

		// 30.01.2018 just for posterior probability computation:
		// double conditional_alignment_total_log_prob = dp_alignment.compute_conditional_alignment_total_log_prob(band_width);
		// myLog << "conditional alignment total LL with gotoh table: " << conditional_alignment_total_log_prob << endl << endl;
	}

	double avg_elapsed_secs_per_dp_alignment = elapsed_secs_in_all_dp_procedures / number_of_computed_dp_alignments;
	myLog << "total time parameter optimization stage (seconds): " << elapsed_secs_param_opt << endl;
	myLog << "total number of computed dp alignments: " << number_of_computed_dp_alignments << endl;
	myLog << "total time (seconds): " << elapsed_secs_in_all_dp_procedures << endl;
	myLog << "average time per dp alignment (seconds): " << avg_elapsed_secs_per_dp_alignment << endl;
	myLog.close();

	// write clean results:
	ofstream myRes;
	myRes.open(result_file_string);

	string how_table_was_chosen;
	string est_table;
	double est_branch_length;
	if (should_gotoh_opt)
	{
		how_table_was_chosen = "gotoh";
		est_table = best_table_gotoh_alignment;
		est_branch_length = gotoh_branch_length_t;
	}
	else
	{
		how_table_was_chosen = "best";
		est_table = best_table_best_alignment;
		est_branch_length = best_branch_length_t;
	}

	if (is_aligned)
	{
		// compute input alignment LL with chosen table:
		vector<vector<size_t>> input_coded_alignment = input_seqs_obj.get_coded_alignment();
		read_chop_tables est_chop_tables_obj(est_table);
		alignment input_alignment(est_chop_tables_obj, input_coded_alignment, est_branch_length, is_jc, quick_jtt_obj);
		double est_table_conditional_input_alignment_log_prob = input_alignment.get_alignment_log_probability_cond_on_anc();

		myRes << "input alignment LL with " << how_table_was_chosen << " table: " << est_table_conditional_input_alignment_log_prob << endl;
		myRes << "input alignment fasta: " << endl;
		print_seqs_fasta(input_seqs, myRes);
	}
	
	myRes << "best alignment " << how_table_was_chosen << " table: " << est_table << endl;
	myRes << "best alignment LL with " << how_table_was_chosen << " table: " << conditional_best_alignment_log_prob << endl;
	myRes << "best alignment fasta: " << endl;
	print_seqs_fasta(best_alignment, myRes);

	for (size_t i = 0; i < (size_t)number_alignments_to_sample; i++)
	{
		vector<string> sampled_alignment;
		double conditional_sampled_alignment_log_prob = 1.0; // according to the table most fitting for the best alignment

		// create chop table object:
		read_chop_tables est_chop_tables_obj(est_table);
		// construct dp object:
		compute_alignment_dp dp_alignment(est_chop_tables_obj, input_seqs_obj, est_branch_length, is_jc, quick_jtt_obj);
		// sample:
		bool should_take_max = false;
		sampled_alignment = dp_alignment.get_sampled_alignment(conditional_sampled_alignment_log_prob, should_take_max, band_width);

		myRes << "sampled alignment " << how_table_was_chosen << " table (according to " << how_table_was_chosen << " alignment): " << est_table << endl;
		myRes << "sampled alignment LL with " << how_table_was_chosen << " table: " << conditional_sampled_alignment_log_prob << endl;
		myRes << "sampled alignment fasta: " << endl;
		print_seqs_fasta(sampled_alignment, myRes);
	}
	
	myRes.close();

	return 0;
}


