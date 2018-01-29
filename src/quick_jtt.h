// code by Tal Pupko and Eli Levy Karin

#ifndef ___QUICK_JTT_H
#define ___QUICK_JTT_H	

#include <cstdio>
#include <cmath>
#include <vector>
#include <random>
#include <iostream>
#include <sstream>
#include <map>
#include <fstream>
#include <string>
using namespace std;

class quick_jtt
{
public:

	quick_jtt(string JTT_coeffs_file)
	{
		read_file_with_cheby_jtt_coef(JTT_coeffs_file);
	} //constructor
	quick_jtt()
	{

	} // empty constructor

	double get_pij_t(const size_t from_aa, const size_t to_aa, const double branch_length_t); // aa's indexing starts at 1 (not 0, like the lib). 0 is always a "-"
	double freq(const size_t aa) const { return _freq[aa - 1]; } // aa's indexing starts at 1 (not 0, like the lib). 0 is always a "-"

private:
	vector<vector<vector<double>>> _chebi_coff;
	vector<double> _freq;
	void read_file_with_cheby_jtt_coef(string JTT_coeffs_file);
};

#endif

#pragma once
