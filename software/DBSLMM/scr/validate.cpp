/*
Deterministic Bayesian Sparse Linear Mixed Model (DBSLMM)
Copyright (C) 2019  Sheng Yang and Xiang Zhou

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <vector>
#include <string>
#include <iostream>
#include <ctime>
#include <math.h>
#include <armadillo>

#include "dtpr.hpp"
#include "validate.hpp"

using namespace std;

VALID::VALID(void) :
	version("0.1"), date("12/08/2019"), year("2019")
{}

void VALID::printHeader(void)
{
	cout << endl;
	cout << "*************************************************************"<< endl;
	cout << "  Deterministic Bayesian Sparse Linear Mixed Model (DBSLMM)  " << endl;
	cout << "  Version " << version << ", " << date << "                  " << endl;
	cout << "  Visit http://www.xzlab.org/software.html For Updates       " << endl;
	cout << "  (C) " << year << " Sheng Yang, Xiang Zhou                  " << endl;
	cout << "  GNU General Public License                                 " << endl;
	cout << "  For Help, Type ./valid -h                                  " << endl;
	cout << "*************************************************************" << endl;
	cout << endl;

	return;
}

void VALID::printHelp(void) {
	cout << " FILE I/O RELATED OPTIONS" << endl;
	cout << " -d        [filename]  " << " specify input the result of DBSLMM." << endl;
	cout << " -s        [filename]  " << " specify input the external summary data." << endl;
	cout << " -r        [filename]  " << " specify input the bfile of reference data." << endl;
	cout << " -mafMax   [num]       " << " specify input the maximium of the difference between reference panel and external data." << endl;
	cout << " -b        [filename]  " << " specify input the block information." << endl;
	cout << " -r2       [num]       " << " specify output r2." << endl;
	return;
}

void VALID::Assign(int argc, char ** argv, PARAM &cPar) {
	
	string str;
	for (int i = 0; i < argc; i++) {

		if (strcmp(argv[i], "--dbslmm") == 0 || strcmp(argv[i], "-d") == 0) {

			if (argv[i + 1] == NULL || argv[i + 1][0] == '-') { continue; }
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.d = str;
		}
		else if (strcmp(argv[i], "--summ") == 0 || strcmp(argv[i], "-s") == 0) {

			if (argv[i + 1] == NULL || argv[i + 1][0] == '-') { continue; }
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.s = str;
		}
		else if (strcmp(argv[i], "--reference") == 0 || strcmp(argv[i], "-r") == 0) {

			if (argv[i + 1] == NULL || argv[i + 1][0] == '-') { continue; }
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.r = str;
		}
		else if (strcmp(argv[i], "--mafMax") == 0 || strcmp(argv[i], "-mafMax") == 0) {

			if (argv[i + 1] == NULL || argv[i + 1][0] == '-') { continue; }
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.mafMax = atof(str.c_str());
		}
		else if (strcmp(argv[i], "--block") == 0 || strcmp(argv[i], "-b") == 0) {

			if (argv[i + 1] == NULL || argv[i + 1][0] == '-') { continue; }
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.b = str;
		}
		else if (strcmp(argv[i], "--R2") == 0 || strcmp(argv[i], "-r2") == 0) {

			if (argv[i + 1] == NULL || argv[i + 1][0] == '-') { continue; }
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.r2 = str;
		}
	}
	return;
}


void VALID::BatchRun(PARAM &cPar) {

	SNPPROC cSP;
	IO cIO;

	// input check
	cout << "Options: " << endl;
	cout << "-d:      " << cPar.d << endl;
	cout << "-s:      " << cPar.s << endl;
	cout << "-r:      " << cPar.r << endl;
	cout << "-mafMax: " << cPar.mafMax << endl;
	cout << "-b:      " << cPar.b << endl;
	cout << "-r2:      " << cPar.r2 << endl;

	// check files
	string ref_fam_str = cPar.r + ".fam"; 
	ifstream dStream(cPar.d.c_str()), sStream(cPar.s.c_str()), rStream(ref_fam_str), bStream(cPar.b.c_str());
	if (cPar.d.size() == 0) {
		cerr << "ERROR: -d is no parameter!" << endl; 
		exit(1); 
	}
	if (cPar.s.size() == 0) {
		cerr << "ERROR: -s is no parameter!" << endl; 
		exit(1); 
	}
	if (cPar.r.size() == 0) {
		cerr << "ERROR: -r is no parameter!" << endl; 
		exit(1); 
	}
	if (cPar.s.size() == 0) {
		cerr << "ERROR: -s is no parameter!" << endl; 
		exit(1); 
	}
	if (!dStream) {
		cerr << "ERROR: " << cPar.d << " dose not exist!" << endl; 
		exit(1);
	}
	if (!sStream) {
		cerr << "ERROR: " << cPar.s << " dose not exist!" << endl; 
		exit(1);
	}
	if (!rStream) {
		cerr << "ERROR: " << cPar.r << " dose not exist!" << endl; 
		exit(1);
	}
	if (!bStream) {
		cerr << "ERROR: " << cPar.b << " dose not exist!" << endl; 
		exit(1);
	}
	
	char separate1[] = " ";
	char separate2[] = "\t";
	// get sample size of reference panel 
	cout << "Reading PLINK FAM file from [" << cPar.r << ".fam]" << endl;
	int n_ref = cIO.getRow(ref_fam_str);
	cout << n_ref << " individuals to be included from reference FAM file." << endl;
	
	// input dbslmm result 
	cout << "Reading DBSLMM result file from [" << cPar.d << "]" << endl;
	vector<SUMMS> dbslmm; 
	cIO.readDBSLMM(cPar.d, separate1, dbslmm); 
	cout << dbslmm.size() << " SNPs in DBSLMM result. " << endl;
	
	// input external summary statistics 
	cout << "Reading summary statistics file from [" << cPar.s << "]" << endl;
	map<string, SUMMS> summ_ext; 
	cIO.readExt(cPar.s, separate1, summ_ext); 
	cout << summ_ext.size() << " SNPs in external result. " << endl;
	
	// match dbslmm and external summary statistics 
	vector<SUMMC> summ_comb; 
	cSP.matchSumm(dbslmm, summ_ext, summ_comb); 
	cout << summ_comb.size() << " SNPs are intersection of DBSLMM and external summary statistics. " << endl;
	clearVector(dbslmm); 
	summ_ext.clear(); 
	
	// input SNP of reference panel
	cout << "Reading PLINK BIM file from [" << cPar.r << ".bim]" << endl; 
	map <string, ALLELEB> ref_bim; 
	bool constr = true; 
	if (abs(cPar.mafMax-1.0) < 1e-10){
		constr = false; 
	}
	cIO.readBim(n_ref, cPar.r, separate2, ref_bim, constr); 
	int num_snp_ref = ref_bim.size(); 
	cout << num_snp_ref << " SNPs to be included from reference BIM file." << endl; 
	
	// match all summ_comb to reference panel 
	vector <SUMMP> summ_comb_r; 
	cSP.matchAll(summ_comb, ref_bim, cPar.mafMax, summ_comb_r); 
	cout << summ_comb_r.size() << " SNPs intersect." << endl;

	// input block file
	vector <BLOCK> block_dat; 
	cIO.readBlock(cPar.b, separate2, block_dat);
	int num_block = block_dat.size(); 
	cout << num_block << " blocks for the chromesome." << endl; 
	
	// calculate for each block 
	vector<int> idv(n_ref);
	for (int i = 0; i < n_ref; i++) idv[i] = 1; 
	string bed_str = cPar.r + ".bed";
	ifstream bed_in(bed_str.c_str(), ios::binary);
	vec nume = zeros<vec>(num_block), deno = zeros<vec>(num_block);
	int count1 = 0, count2 = 0; 
	for (int b = 0; b < num_block; ++b){
		int num_snp_b = 0; 
		long start = block_dat[b].start;
		long end = block_dat[b].end;
		for (size_t j = count1; j < summ_comb_r.size(); j++){
			if (summ_comb_r[j].ps >= start && summ_comb_r[j].ps < end) { 
				num_snp_b++; 
				count1++; 
			}else{
				break;
			}
		}
		vec z1_b = zeros<vec>(num_snp_b), z2_b = zeros<vec>(num_snp_b); 
		mat geno_b = zeros<mat>(n_ref, num_snp_b); 
		int count = 0; 
		for (size_t j = count2; j < summ_comb_r.size(); j++){
			if (summ_comb_r[j].ps >= start && summ_comb_r[j].ps < end) { 
				z1_b(count) = summ_comb_r[j].z1; 
				z2_b(count) = summ_comb_r[j].z2; 
				vec geno = zeros<vec>(n_ref);
				double maf = 0.0; 
				cIO.readSNPIm(summ_comb_r[j].pos, n_ref, idv, bed_in, geno, maf);
				cSP.nomalizeVec(geno);
				geno_b.col(count) = geno;
				count2++; 
				count++; 
			}else{
				break;
			}
		}
		mat SIGMA_b = geno_b.t() * geno_b; 
		SIGMA_b /= (double)n_ref; 
		nume[b] = as_scalar(z1_b.t() * z2_b); 
		deno[b] = as_scalar(z1_b.t() * SIGMA_b * z1_b); 
	}
	string r2_str = cPar.r2 + ".txt"; 
	ofstream r2Fout(r2_str.c_str());
	for (int i = 0; i < num_block; ++i) 
		r2Fout << nume[i] << " " << deno[i] << endl; 
	
	return; 
}