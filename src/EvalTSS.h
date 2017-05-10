/*************************************************************************
EvalTSS.h
pfm, 2015, 2016

>>> SOURCE LICENSE >>>
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation (www.fsf.org); either version 2 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

A copy of the GNU General Public License is available at
http://www.fsf.org/licensing/licenses
>>> END OF LICENSE >>>
*************************************************************************/
#ifndef EVALTSS_H
#define EVALTSS_H
#include <iostream>
#include <fstream>
#include <sstream>
#include <cassert>
#include <math.h>
#include <cfloat>
#include <map>
#include <string>
#include "TSS.h"

#define NL 0
#define INFTY 1E10

using namespace std;

class EvalTSS {
public:
	std::vector<std::vector<int>> confusionMatrix;
	double precision, recall, sensitivity, specificity, accuracy, F1, errorRate, globalErrorRate;
	int TP, TN, FP, FN, nFN, RPT, Lref, N;
	double alignmentScore, relativeAlignmentScore;
	double relativeLatency, relativeRPT;
	double meanDuration;
	unsigned int nbLabs;

	EvalTSS();
	EvalTSS(int N);

	void add_cm(std::vector<std::vector<int>> cm);
	static int  countDistinctSymbols(std::vector<int> ls);
	void eval(TSS ref, TSS pred);
	void eval(std::string &fileRef, std::string &filePred);
	void evalDP(std::string &fileRef, std::string &filePred);

	void evalIsolatedPatterns(std::string &fileRef, std::string &filePred);
	void saveFile(std::string file, bool append, int nTEST);
	void saveFileA(std::string fname, bool append, int nTEST);

	void saveStaticLabels(std::string filename, std::vector<int> labs);

	void readFileConfMat(std::string fname);
	void saveFileConfMat(std::string fname);
	std::vector<int> loadLabels(std::string &filename, int skip);
    	void SeqAlignWithBacktrace(TSS tr, TSS tp);

	void print();
	void print(int beg);
	void evaluateMetricsFromConfusionMatrix();
	void evaluateMetricsFromConfusionMatrix(int beg);
	void evaluateMetricsFromConfusionMatrix_1vsAll(int beg, int ilab);
	void initializeConfusionMatrix(int N);
    };

#endif
