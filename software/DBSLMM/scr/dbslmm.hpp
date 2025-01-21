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

#ifndef __DBSLMM_H__
#define __DBSLMM_H__

#include <vector>
#include <string>
#include <iostream>

#include "dtpr.hpp"

class PARAM {
public:
	//parameters
	string s; 
	string l;
	string r;
	int n;
	double mafMax; 
	int nsnp;
	string b;
	double h;
	int t;
	string eff;
};

class DBSLMM {
public:
	//parameters
	string version;
	string date;
	string year;

	//constructor
	DBSLMM(void);

	//functions
	void printHeader(void);
	void printHelp(void);
	void Assign(int argc, char ** argv, PARAM &cPar);
	void BatchRun(PARAM &cPar);
};


#endif