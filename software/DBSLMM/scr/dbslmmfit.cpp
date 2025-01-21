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

#define ARMA_DONT_USE_WRAPPER

#include <iostream>
#include <vector>
#include <armadillo>
#include <string>
#include "omp.h"

#include "dtpr.hpp"
#include "dbslmmfit.hpp"

using namespace std;
using namespace arma;

// estimate large and small effect
int DBSLMMFIT::est(int n_ref, int n_obs, double sigma_s, int num_block, vector<int> idv, string bed_str,
				 vector <INFO> info_s, vector <INFO> info_l, int thread, 
                 vector <EFF> &eff_s, vector <EFF> &eff_l){
	
	// get the maximum number of each block
	int count_s = 0, count_l = 0;
	vec num_s = zeros<vec>(num_block), num_l = zeros<vec>(num_block); 
	for (int i = 0; i < num_block; i++) {
		for (size_t j = count_s; j < info_s.size(); j++) {
			if(info_s[j].block == i){ 
				num_s(i) += 1; 
				count_s++;
			}else{
				break;
			}
		}
		for (size_t j = count_l; j < info_l.size(); j++) {
			if(info_l[j].block == i){ 
				num_l(i) += 1; 
				count_l++;
			}else{
				break;
			}
		}
	}
	count_l = count_s = 0; // reset
	
	double len_l = num_l.max(); 
	double len_s = num_s.max(); 
	
	int B = 0;
	int B_MAX = 60;
	if (num_block < 60){
		B_MAX = num_block; 
	}

	// pseudo INFO
	INFO info_pseudo; 
	info_pseudo.snp = "rs"; 
	info_pseudo.ps = 0; 
	info_pseudo.pos = 0; 
	info_pseudo.block = 0; 
	info_pseudo.a1 = "Y"; 
	info_pseudo.maf = 0; 
	info_pseudo.z = 0; 
	info_pseudo.P = 0; 
	
	// loop 
	// vector < vector <INFO*> > info_s_Block(B_MAX, vector <INFO*> ((int)len_s)), info_l_Block(B_MAX, vector <INFO*> ((int)len_l));
	vector < vector <INFO> > info_s_Block(B_MAX, vector <INFO> ((int)len_s)), info_l_Block(B_MAX, vector <INFO> ((int)len_l));
	vector < vector <EFF> > eff_s_Block(B_MAX, vector <EFF> ((int)len_s)), eff_l_Block(B_MAX, vector <EFF> ((int)len_l));
	vector <int> num_s_vec, num_l_vec;
	for (int i = 0; i < num_block; ++i) {
		// small effect SNP information
		vector <INFO> info_s_block; 
		for (size_t j = count_s; j < info_s.size(); j++) {
			if(info_s[j].block == i){ 
				info_s_block.push_back(info_s[j]);
				count_s++;
			}else{
				break;
			}
		}
		for (size_t k = 0; k < info_s_block.size(); k++)
			// info_s_Block[B][k] = &info_s_block[k]; 
			info_s_Block[B][k] = info_s_block[k]; 
		num_s_vec.push_back((int)num_s(i));
		
		// large effect SNP information
		if (num_l(i) == 0){
			// info_l_Block[B][0] = &info_pseudo;
			info_l_Block[B][0] = info_pseudo;
		}else{
			vector <INFO> info_l_block;
			for (size_t j = count_l; j < info_l.size(); j++) {
				if(info_l[j].block == i){ 
					info_l_block.push_back(info_l[j]);
					count_l++;
				}else{
					break;
				}
			}
			for (size_t k = 0; k < info_l_block.size(); k++)
				// info_l_Block[B][k] = &info_l_block[k]; 
				info_l_Block[B][k] = info_l_block[k]; 
		}
		num_l_vec.push_back((int)num_l(i));

		B++;
		if (B == B_MAX || i + 1 == num_block) { // process the block of SNPs using multi-threading
			
			omp_set_num_threads(thread);
#pragma omp parallel for schedule(dynamic)
			for (int b = 0; b < B; b++){
				calcBlock(n_ref, n_obs, sigma_s, idv, bed_str, info_s_Block[b], info_l_Block[b],
						  num_s_vec[b], num_l_vec[b], eff_s_Block[b], eff_l_Block[b]);
			}
			// eff of small effect SNPs
			for (int r = 0; r < B; r++) {
				for (int l = 0; l < num_s_vec[r]; l++){
					eff_s.push_back(eff_s_Block[r][l]);
				}
			}
			// eff of large effect SNPs
			for (int r = 0; r < B; r++){
				for (int l = 0; l < num_l_vec[r]; l++){
					eff_l.push_back(eff_l_Block[r][l]);
				}
			}
			B = 0;
			num_l_vec.clear(); 
			num_s_vec.clear();
		}
	}
	return 0;
}

// estimate only small effect
int DBSLMMFIT::est(int n_ref, int n_obs, double sigma_s, int num_block, vector<int> idv, string bed_str,
				 vector <INFO> info_s, int thread, 
				 vector <EFF> &eff_s){
	
	// get the maximum number of each block
	int count_s = 0;
	vec num_s = zeros<vec>(num_block); 
	for (int i = 0; i < num_block; i++) {
		for (size_t j = count_s; j < info_s.size(); j++) {
			if(info_s[j].block == i){ 
				num_s(i) += 1; 
				count_s++;
			}else{
				break;
			}
		}
	}
	count_s = 0; // reset
	
	double len_s = num_s.max(); 
	
	int B = 0;
	int B_MAX = 60;
	if (num_block < 60){
		B_MAX = num_block; 
	}

	// pseudo INFO
	INFO info_pseudo; 
	info_pseudo.snp = "rs"; 
	info_pseudo.ps = 0; 
	info_pseudo.pos = 0; 
	info_pseudo.block = 0; 
	info_pseudo.a1 = "Y"; 
	info_pseudo.maf = 0; 
	info_pseudo.z = 0; 
	info_pseudo.P = 0; 
	
	// loop 
	// vector < vector <INFO*> > info_s_Block(B_MAX, vector <INFO*> ((int)len_s));
	vector < vector <INFO> > info_s_Block(B_MAX, vector <INFO> ((int)len_s));
	vector < vector <EFF> > eff_s_Block(B_MAX, vector <EFF> ((int)len_s));
	vector <int> num_s_vec;
	for (int i = 0; i < num_block; ++i) {
		// small effect SNP information
		vector <INFO> info_s_block; 
		for (size_t j = count_s; j < info_s.size(); j++) {
			if(info_s[j].block == i){ 
				info_s_block.push_back(info_s[j]);
				count_s++;
			}else{
				break;
			}
		}
		for (size_t k = 0; k < info_s_block.size(); k++)
			// info_s_Block[B][k] = &info_s_block[k]; 
			info_s_Block[B][k] = info_s_block[k]; 
		num_s_vec.push_back((int)num_s(i));
		
		B++;
		if (B == B_MAX || i + 1 == num_block) { // process the block of SNPs using multi-threading
			
			omp_set_num_threads(thread);
#pragma omp parallel for schedule(dynamic)
			for (int b = 0; b < B; b++){
				calcBlock(n_ref, n_obs, sigma_s, idv, bed_str, info_s_Block[b],
						  num_s_vec[b], eff_s_Block[b]);
			}
			// eff of small effect SNPs
			for (int r = 0; r < B; r++) {
				for (int l = 0; l < num_s_vec[r]; l++){
					eff_s.push_back(eff_s_Block[r][l]);
				}
			}
			B = 0;
			num_s_vec.clear();
		}
	}
	return 0;
}

// estimate large and small effect for each block
int DBSLMMFIT::calcBlock(int n_ref, int n_obs, double sigma_s, vector<int> idv, string bed_str, 
						vector <INFO> info_s_block_full, vector <INFO> info_l_block_full, 
						int num_s_block, int num_l_block, 
						vector <EFF> &eff_s_block, vector <EFF> &eff_l_block){
	SNPPROC cSP;
	IO cIO; 
	ifstream bed_in(bed_str.c_str(), ios::binary);
	
	// INFO small effect SNPs 
	// vector <INFO*> info_s_block(num_s_block);
	// for (int i = 0; i < num_s_block; i++)
		// info_s_block[i] = info_s_block_full[i];
	vector <INFO> info_s_block(num_s_block);
	for (int i = 0; i < num_s_block; i++)
		info_s_block[i] = info_s_block_full[i];
	// z_s
	vec z_s = zeros<vec>(num_s_block); 
	for (int i = 0; i < num_s_block; i++) 
		// z_s(i) = info_s_block[i]->z;
		z_s(i) = info_s_block[i].z;
	// small effect genotype matrix
	mat geno_s = zeros<mat>(n_ref, num_s_block);
	for (int i = 0; i < num_s_block; ++i) {
		vec geno = zeros<vec>(n_ref);
		double maf = 0.0; 
		// cIO.readSNPIm(info_s_block[i]->pos, n_ref, idv, bed_in, geno, maf);
		cIO.readSNPIm(info_s_block[i].pos, n_ref, idv, bed_in, geno, maf);
		cSP.nomalizeVec(geno);
		geno_s.col(i) = geno;
	}
	// pseudo EFF
	EFF eff_pseudo; 
	eff_pseudo.snp = "rs"; 
	eff_pseudo.a1 = "Y"; 
	eff_pseudo.maf = 0.0; 
	eff_pseudo.beta = 0.0; 
	
	vec beta_s = zeros<vec>(num_s_block); 
	// INFO large effect SNPs 
	if (num_l_block != 0){
		// vector <INFO*> info_l_block(num_l_block);
		vector <INFO> info_l_block(num_l_block);
		for (int i = 0; i < num_l_block; i++) 
			info_l_block[i] = info_l_block_full[i];
		// z_l
		vec z_l = zeros<vec>(num_l_block); 
		for(int i = 0; i < num_l_block; i++) 
			// z_l(i) = info_l_block[i]->z;
			z_l(i) = info_l_block[i].z;


		// large effect matrix
		mat geno_l = zeros<mat>(n_ref, num_l_block);
		for (int i = 0; i < num_l_block; ++i) {
			vec geno = zeros<vec>(n_ref);
			double maf = 0.0; 
			// cIO.readSNPIm(info_l_block[i]->pos, n_ref, idv, bed_in, geno, maf);
			cIO.readSNPIm(info_l_block[i].pos, n_ref, idv, bed_in, geno, maf);
			cSP.nomalizeVec(geno);
			geno_l.col(i) = geno;
		}
		// estimation
		vec beta_l = zeros<vec>(num_l_block); 
		estBlock(n_ref, n_obs, sigma_s, geno_s, geno_l, z_s, z_l, beta_s, beta_l);
		// summary 
		for(int i = 0; i < num_l_block; i++) {
			EFF eff_l; 
			// eff_l.snp = info_l_block[i]->snp;
			// eff_l.a1 = info_l_block[i]->a1;
			// eff_l.maf = info_l_block[i]->maf;
			eff_l.snp = info_l_block[i].snp;
			eff_l.a1 = info_l_block[i].a1;
			eff_l.maf = info_l_block[i].maf;
			eff_l.beta = beta_l(i);
			eff_l_block[i] = eff_l;
		}
	}
	else{
		// estimation
		estBlock(n_ref, n_obs, sigma_s, geno_s, z_s, beta_s);
		eff_l_block[0].snp = eff_pseudo.snp;
		eff_l_block[0].a1 = eff_pseudo.a1;
		eff_l_block[0].maf = eff_pseudo.maf;
		eff_l_block[0].beta = eff_pseudo.beta;
	}
	// output small effect
	for(int i = 0; i < num_s_block; i++) {
		EFF eff_s; 
		// eff_s.snp = info_s_block[i]->snp;
		// eff_s.a1 = info_s_block[i]->a1;
		// eff_s.maf = info_s_block[i]->maf;
		eff_s.snp = info_s_block[i].snp;
		eff_s.a1 = info_s_block[i].a1;
		eff_s.maf = info_s_block[i].maf;
		eff_s.beta = beta_s(i); 
		eff_s_block[i] = eff_s;
	}
	return 0; 
}

// estimate only small effect for each block
int DBSLMMFIT::calcBlock(int n_ref, int n_obs, double sigma_s, vector<int> idv, string bed_str, 
						vector <INFO> info_s_block_full, int num_s_block, 
						vector <EFF> &eff_s_block){
	SNPPROC cSP;
	IO cIO; 
	ifstream bed_in(bed_str.c_str(), ios::binary);
	
	// INFO small effect SNPs 
	// vector <INFO*> info_s_block(num_s_block);
	// for (int i = 0; i < num_s_block; i++)
		// info_s_block[i] = info_s_block_full[i];
	vector <INFO> info_s_block(num_s_block);
	for (int i = 0; i < num_s_block; i++)
		info_s_block[i] = info_s_block_full[i];
	// z_s
	vec z_s = zeros<vec>(num_s_block); 
	for (int i = 0; i < num_s_block; i++) 
		// z_s(i) = info_s_block[i]->z;
		z_s(i) = info_s_block[i].z;
	// small effect genotype matrix
	mat geno_s = zeros<mat>(n_ref, num_s_block);
	for (int i = 0; i < num_s_block; ++i) {
		vec geno = zeros<vec>(n_ref);
		double maf = 0.0; 
		// cIO.readSNPIm(info_s_block[i]->pos, n_ref, idv, bed_in, geno, maf);
		cIO.readSNPIm(info_s_block[i].pos, n_ref, idv, bed_in, geno, maf);
		cSP.nomalizeVec(geno);
		geno_s.col(i) = geno;
	}
	
	// estimation
	vec beta_s = zeros<vec>(num_s_block); 
	estBlock(n_ref, n_obs, sigma_s, geno_s, z_s, beta_s);
	
	// output small effect
	for(int i = 0; i < num_s_block; i++) {
		EFF eff_s; 
		// eff_s.snp = info_s_block[i]->snp;
		// eff_s.a1 = info_s_block[i]->a1;
		// eff_s.maf = info_s_block[i]->maf;
		eff_s.snp = info_s_block[i].snp;
		eff_s.a1 = info_s_block[i].a1;
		eff_s.maf = info_s_block[i].maf;
		eff_s.beta = beta_s(i); 
		eff_s_block[i] = eff_s;
	}
	return 0; 
}

// solve the equation Ax=b, x is a variables
vec DBSLMMFIT::PCGv(mat A, vec b, size_t maxiter, const double tol){
	vec dA = A.diag();
	// stable checking, using available func to speed up
	for(size_t i = 0; i< dA.n_elem; i++){
		if(dA[i] == 0)
			dA[i] = 1e-4;
	}// end for i
	vec Minv = 1.0/dA;
	// initialize
	vec x = zeros<vec>(b.n_elem);
	vec r = zeros<vec>(b.n_elem);
	vec r1 = zeros<vec>(b.n_elem);
	vec z1 = zeros<vec>(b.n_elem);
	r = b;
	vec z = Minv % r;
	vec p = z;
	size_t iter = 0;
	double sumr2 = norm(r, 2);
	// PCG main loop 
	while( sumr2 > tol && iter < maxiter){
		iter += 1;
		// move direction
		vec Ap = A*p;
		// step size
		double a = dot(r,z)/dot(p,Ap);
		// move
		x = x + a * p;
		r1 = r - a * Ap;
		z1 = Minv % r1;
		double bet = dot(z1, r1)/dot(z, r);
		p = z1 + bet * p;
		z = z1;
		r = r1;
		sumr2 = norm(r, 2);
	}// end while loop
	if (iter >= maxiter){
		cerr << "ERROR: Matrix is Singular!" << endl;
	}
	return(x);
}// end function

mat DBSLMMFIT::PCGm(mat A, mat B, size_t maxiter, const double tol){
	
	size_t n_iter = B.n_cols;
	mat x = zeros<mat>(A.n_rows, n_iter);
	for (size_t i = 0; i < n_iter; i++){
		x.col(i) = PCGv(A, B.col(i), maxiter, tol);
	}// end for loop
	return(x);
}// end function

int DBSLMMFIT::estBlock(int n_ref, int n_obs, double sigma_s, mat geno_s, mat geno_l, vec z_s, vec z_l, vec &beta_s, vec &beta_l) {
	
	// LD matrix 
	// mat SIGMA_ls = geno_l.t() * geno_s; 
	// SIGMA_ls /= (double)n_ref; 
	// mat SIGMA_ll = geno_l.t() * geno_l; 
	// SIGMA_ll /= (double)n_ref;
	// mat SIGMA_ss = geno_s.t() * geno_s; 
	// SIGMA_ss /= (double)n_ref;
	
	double tau = 0.8;
	mat SIGMA_ls = geno_l.t() * geno_s;
	SIGMA_ls *= tau/(double)n_ref;
	mat SIGMA_ll = geno_l.t() * geno_l; 
	SIGMA_ll *= tau/(double)n_ref;
	vec DIAG_l(geno_l.n_cols, fill::ones);
	DIAG_l *= (1.0-tau);
	SIGMA_ll += diagmat(DIAG_l);
	mat SIGMA_ss = geno_s.t() * geno_s; 
	SIGMA_ss *= tau/(double)n_ref;
	vec DIAG_s(geno_s.n_cols, fill::ones);
	DIAG_s *= (1.0-tau);
	SIGMA_ss += diagmat(DIAG_s);

	// beta_l
	SIGMA_ss.diag() += 1.0 / (sigma_s * (double)n_obs);
	mat SIGMA_ss_inv_SIGMA_sl = PCGm(SIGMA_ss, SIGMA_ls.t(), 1000, 1e-7);
	mat SIGMA_ls_SIGMA_ss_inv_SIGMA_sl = - SIGMA_ls * SIGMA_ss_inv_SIGMA_sl;
	SIGMA_ls_SIGMA_ss_inv_SIGMA_sl += SIGMA_ll;
	vec SIGMA_ss_inv_z_s = PCGv(SIGMA_ss, z_s, 1000, 1e-7);
	vec SIGMA_ls_SIGMA_ss_inv_z_s = -SIGMA_ls * SIGMA_ss_inv_z_s;
	SIGMA_ls_SIGMA_ss_inv_z_s += z_l;
	beta_l = PCGv(SIGMA_ls_SIGMA_ss_inv_SIGMA_sl, SIGMA_ls_SIGMA_ss_inv_z_s, 1000, 1e-7);
	beta_l /= sqrt(n_obs);

	// beta_s
	SIGMA_ss_inv_z_s *= sqrt(n_obs);
	vec SIGMA_ss_inv_SIGMA_sl_beta_l = (double)n_obs * SIGMA_ss_inv_SIGMA_sl * beta_l; 
	vec SIGMA_ss_inv_z_s_SIGMA_sl_beta_l = SIGMA_ss_inv_z_s - SIGMA_ss_inv_SIGMA_sl_beta_l; 
	SIGMA_ss.diag() -= 1.0 / (sigma_s * (double)n_obs);
	vec SIGMA_ss_z_s_SIGMA_sl_beta_l = SIGMA_ss * SIGMA_ss_inv_z_s_SIGMA_sl_beta_l; 
	beta_s = sqrt(n_obs) * z_s - (double)n_obs * SIGMA_ls.t() * beta_l - SIGMA_ss_z_s_SIGMA_sl_beta_l; 
	beta_s *= sigma_s;
	
	return 0; 
}

int DBSLMMFIT::estBlock(int n_ref, int n_obs, double sigma_s, mat geno_s, vec z_s, vec &beta_s) {
	
	// LD matrix 
	// mat SIGMA_ss = geno_s.t() * geno_s; 
	// SIGMA_ss /= (double)n_ref;

	double tau = 0.8;
	mat SIGMA_ss = geno_s.t() * geno_s; 
	SIGMA_ss *= tau/(double)n_ref;
	vec DIAG_s(geno_s.n_cols, fill::ones);
	DIAG_s *= (1.0-tau);
	SIGMA_ss += diagmat(DIAG_s);
	
	// beta_s
	SIGMA_ss.diag() += 1.0 / (sigma_s * (double)n_obs);
	vec SIGMA_ss_inv_z_s = PCGv(SIGMA_ss, z_s, 1000, 1e-7);
	SIGMA_ss.diag() -= 1.0 / (sigma_s * (double)n_obs);
	vec SIGMA_ss_SIGMA_ss_inv_z_s = SIGMA_ss * SIGMA_ss_inv_z_s;
	vec z_s_SIGMA_ss_SIGMA_ss_inv_SIGMA_sl = z_s - SIGMA_ss_SIGMA_ss_inv_z_s; 
	beta_s = sqrt(n_obs) * sigma_s * z_s_SIGMA_ss_SIGMA_ss_inv_SIGMA_sl; 
	
	return 0; 
}
