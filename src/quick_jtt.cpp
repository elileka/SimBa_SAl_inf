// code by Tal Pupko and Eli Levy Karin

#include "quick_jtt.h"

double quick_jtt::get_pij_t(const size_t from_aa, const size_t to_aa, const double branch_length_t)
{
	double d = 0.0, dd = 0.0, sv, y, y2, pij_t_prob;
	double _leftRange = 2.0;
	double _rightRange = 0.0;
	int _usingNumberOfCoef = 13;
	int j;

	size_t from_aa_ind = from_aa - 1; // aa's indexing starts at 1 (not 0, like the lib). 0 is always a "-"
	size_t to_aa_ind = to_aa - 1; // aa's indexing starts at 1 (not 0, like the lib). 0 is always a "-"

	y2 = 2.0*(y = (2.0*branch_length_t - _leftRange - _rightRange) / (_rightRange - _leftRange));
	for (j = _usingNumberOfCoef; j>0; j--) 
	{
		sv = d;
		d = y2*d - dd + _chebi_coff[from_aa_ind][to_aa_ind][j];
		dd = sv;
	}
	pij_t_prob = y*d - dd + 0.5*_chebi_coff[from_aa_ind][to_aa_ind][0];
	
	return pij_t_prob;
}

void quick_jtt::read_file_with_cheby_jtt_coef(string JTT_coeffs_file)
{
	fstream f(JTT_coeffs_file.c_str());
	_chebi_coff.resize(20);
	for (size_t i = 0; i < _chebi_coff.size(); ++i) {
		_chebi_coff[i].resize(20);
	}
	for (size_t i = 0; i < _chebi_coff.size(); ++i) {
		for (size_t j = 0; j < _chebi_coff[i].size(); ++j) {
			_chebi_coff[i][j].resize(60);
		}
	}
	for (size_t i = 0; i < _chebi_coff.size(); ++i) {
		for (size_t j = 0; j < _chebi_coff[i].size(); ++j) {
			for (size_t k = 0; k < _chebi_coff[i][j].size(); ++k) {
				f >> _chebi_coff[i][j][k];
			}
		}
	}

	_freq.resize(20);
	for (size_t i = 0; i < 20; ++i) {
		f >> _freq[i];
	}
	f.close();
}






/*
// Tal produced the table of coefficients based on the following code:


#include "pijAccelerator.h"
#include "chebyshevAccelerator.h"
#include "readDatMatrix.h"
int generate_a_file_with_cheby_jtt_coef() {

// GENERATING OF A FILE WITH CHEBY COEFFICIENT BASED ON JTT
// IT USES THE LIB FOR DOING THAT
// HERE IT IS FOR THE JTT MATRIX
// FOR distances (t) between 0 and 2
// THE LAST 20 NUMBERS ARE THE FREQUENCIES

replacementModel* rpc = new pupAll(datMatrixHolder::jones);
//pijAccelerator *pijAcc = new chebyshevAccelerator(rpc);
chebyshevAccelerator ca(rpc,20,60,13,0,2.0);
VVVdouble chebi_coff = ca.chebi_coff;

ofstream of("chebi_coef.txt");
for (size_t i = 0; i < chebi_coff.size(); ++i) {
for (size_t j = 0; j < chebi_coff[i].size(); ++j) {
for (size_t k = 0; k < chebi_coff[i][j].size(); ++k) {
of << setprecision(40)<< chebi_coff[i][j][k] << " ";
}
}
}

for (size_t i = 0; i < 20; ++i) {
of << setprecision(40) << ca.freq(i) << " ";
}

of.close();
return 0;
}



int validateNewPijImplementationCompareToLib() {
generate_a_file_with_cheby_jtt_coef();
replacementModel* rpc = new pupAll(datMatrixHolder::jones);
chebyshevAccelerator ca(rpc, 20, 60, 13, 0, 2.0);

cout << " ------------- new pij implementation ---------------" << endl;
QUICK_PIJ pij_t("chebi_coef.txt");

int fromAA = 5;
int toAA = 3;
double dist = 0.001;;
cout << "new chebyshev pijt = " << pij_t.Pij_t(fromAA, toAA, dist) << endl;
cout << "new chebyshev freq[5] = " << pij_t.freq(5) << endl;
cout << " ------------- old pij implementation ---------------" << endl;
pijAccelerator *pijAcc = new chebyshevAccelerator(rpc);
//cout << "_totalNumOfCoef = " << ca._totalNumOfCoef << endl;;
//cout << "_usingNumberOfCoef = " << ca._usingNumberOfCoef << endl;
//cout << "_rightRange = " << ca._rightRange << endl;
//cout << "_leftRange = " << ca._leftRange << endl;;

cout << "old chebyshev pijt = " << ca.Pij_t(fromAA, toAA, dist) << endl;
cout << "old chebyshev freq[5] = " << ca.freq(5) << endl;

cout << " ------------- regular pij implementation ---------------" << endl;
cout << "regular pijt = " << rpc->Pij_t(fromAA, toAA, dist) << endl;
cout << "regular freq[5] = " << rpc->freq(5) << endl;

return 0;
}

*/