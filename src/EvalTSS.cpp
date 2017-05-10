/*************************************************************************
EvalTSS.cpp
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
#include <memory>
#include <iostream>
#include <fstream>
#include <sstream>
#include "EvalTSS.h"

using namespace std;

double __max(double a, double b){
	if(a>b)
		return a;
	else
		return b;
}
double __min(double a, double b){
	if(a>b)
		return b;
	else
		return a;
}


void EvalTSS::add_cm(std::vector<std::vector<int>> cm){
	for(unsigned int i=0; i<cm.size(); i++){
		for(unsigned int j=0; j< cm[i].size(); j++){
			confusionMatrix[i][j]+=cm[i][j];
		}
	}
}

int  EvalTSS::countDistinctSymbols(std::vector<int> ls){
	std::map<int,int> tab;
	unsigned int LL=ls.size();
	int l=0;
	for (unsigned int i=0; i<LL; i++){
		if(tab.find(ls[i]) == tab.end()){
			tab[(int)ls[i]]=l;
			l++;
		}
	}
	return l;
}




void EvalTSS::SeqAlignWithBacktrace(TSS tr, TSS tp) {
//    Use dynamic programming to find a min-cost path through matrix D.
	int r = tr.data.size();
	int c = tp.data.size();
	double dist;
	std::map<int,int> mapLabels_r=tr.getLabels();
	std::map<int,int> mapLabels_p=tp.getLabels();
	nbLabs=__max(mapLabels_r.size(), mapLabels_p.size());

	initializeConfusionMatrix(nbLabs);

// costs and back-trace matrices 
	double **D = (double **) calloc(r + 1, sizeof(double*));
	int **phi = (int **) calloc(r + 1, sizeof(int*));
	for (int i = 0; i <= r; i++){
		D[i] = (double *) calloc(c + 1, sizeof(double));
		phi[i] = (int *) calloc(c + 1, sizeof(int));
	}

	for (int i = 1; i <= r; i++)
		for (int j = 1; j <= c; j++) {
			dist=0;
			if((tr.data[i-1].beg>tp.data[j - 1].end) || (tr.data[i-1].end<tp.data[j - 1].beg)) //no intersection
					dist = INFTY; // forbid matching
			else if (tr.data[i - 1].lab == tp.data[j - 1].lab){ //intersection and match
				dist=1.0-(__min(tr.data[i-1].end, tp.data[j-1].end)-__max(tr.data[i-1].beg, tp.data[j-1].beg))/
				(__max(tr.data[i-1].end, tp.data[j-1].end)-__min(tr.data[i-1].beg, tp.data[j-1].beg));
				}
			else if(tp.data[j - 1].lab!=NL && tr.data[j - 1].lab!=NL) //intersection and mismatch
				dist = 2;
			else dist=INFTY;

			D[i][j] = dist;
		}
	for (int i = 1; i <= r; i++){
		D[i][0] = 2*i;
		phi[i][0]=2;
	}
	for (int j = 1; j <= c; j++){
		D[0][j] = 2*j;
		phi[0][j]=3;
	}
	D[0][0] = 0.0;


	double dmax;
	int tb;

	for (int i = 1; i <= r; i++) {
		for (int j = 1; j <= c; j++) {
			dmax = D[i - 1][j - 1] + D[i][j];
			tb = 1;
			if (dmax > D[i - 1][j] + 2) {
				dmax = D[i - 1][j] + 2;
				tb = 2;
			}
			if (dmax > D[i][j - 1] + 2) {
				dmax = D[i][j - 1] + 2;
				tb = 3;
			}
			D[i][j] = dmax;
			phi[i][j] = tb;
		}
	}


// Traceback

	int ii = r, i;
	int jj = c, j;

	cout << "GROUND TRUTH" << "\t\t" << "PREDICTION" << "\t\t EVENTS" << endl;
    bool match=false;
    bool mismatch=false;
    bool repet=false;
    alignmentScore=0;
    RPT=0;

    int nbMatch=0;
	while ((ii > 0) || (jj > 0)) {
           i=ii; if(ii-1<0) i=1;
	   j=jj; if(jj-1<0) j=1;
	   cout << i-1 << " (" << tr.data[i-1].lab << "|" << tr.data[i-1].beg << "-" << tr.data[i-1].end << ")\t\t";
	   cout << j-1 << " (" << tp.data[j-1].lab << "|" << tp.data[j-1].beg << "-" << tp.data[j-1].end << ")\t" ;
	   tb = phi[ii][jj];
           //cout << D[i][j] << ":" << tb << ":";
		if (tb == 1) { //match
			if (tr.data[i-1].lab==tp.data[j-1].lab){
				if(tr.data[i-1].lab!=NL){
					alignmentScore+=1.0-(__min(tr.data[i-1].end, tp.data[j-1].end)-__max(tr.data[i-1].beg, tp.data[j-1].beg))/
									(__max(tr.data[i-1].end, tp.data[j-1].end)-__min(tr.data[i-1].beg, tp.data[j-1].beg));
					cout << " correct match type-1\n";
					nbMatch++;
					relativeLatency+=-(tr.data[i-1].end+tr.data[i-1].beg)/2.0+(tp.data[j-1].end+tp.data[j-1].beg)/2.0;
					meanDuration+=tp.data[j-1].end-tp.data[j-1].beg;
					match=true;
				}
				else{
					cout << " NL match type-1\n";
				}
				//confusionMatrix[tp.data[j-1].lab][tr.data[i-1].lab]+=1;
			}
			else {
				alignmentScore+=1;
				if(tp.data[j-1].lab!=NL&&tr.data[i-1].lab!=NL){
					mismatch=true;
					cout << " mismatch (1 FP and 1 FN) type-1\n";
				}
				else if (tr.data[j-1].lab==0){
					cout << " FP 1\n";
				}
				else{
					cout << " FN 1\n";
				}
				//confusionMatrix[tp.data[j-1].lab][tr.data[i-1].lab]+=1;
			}
			confusionMatrix[tp.data[j-1].lab][tr.data[i-1].lab]+=1;
			ii = ii - 1;
			jj = jj - 1;
		} // match
		else if (tb == 2) {//suppress in GT
			if(tr.data[i-1].lab!=NL){
				alignmentScore+=1;
				cout << " suppress in GT (FN) type-2\n";
				confusionMatrix[0][tr.data[i-1].lab]+=1; //FN = mismatch with class 0
			}
			else
				cout << " suppress NL in GT type-2\n";
			ii = ii - 1;
			match=false;
		} // suppress in GT
		else if (tb == 3) {//suppress in pred
			if (tp.data[j-1].lab !=NL){
				alignmentScore+=1;
				cout << " suppress in pred (FP) type-3\n";
				mismatch=true;
				confusionMatrix[tp.data[j-1].lab][0]+=1;
			}
			else{
				cout << " suppress NL in pred type-3 " << endl;
				
			}
			jj = jj - 1;
		} else
			cout << "error !!!!!!!!!!!!!!!!!!!!!	error\n";
	} //while
	relativeAlignmentScore=(float)alignmentScore/(float)r;
	relativeLatency/=nbMatch;
	meanDuration/=nbMatch;
	for (int i = 0; i < r; i++) {
		free(phi[i]);
		free(D[i]);
	}
	free(D[r]);
	free(D);
	free(phi);
} //SeqAlignWithBacktrace


void EvalTSS::print(){
	print(0);
}
void EvalTSS::print(int beg){


    //cout << "N=" << N << " TP=" << TP <<" FP=" << FP <<" TN=" << TN << " FN=" << FN << endl;
	cout << "Macro Error Rate (N-sum_i(FPi+FNi))/N= " << 100*globalErrorRate << " (%)" << endl;
	cout << "Micro Average Precision avg(TPi/(TPi+FPi))= " << 100*precision << " (%)" << endl;
	cout << "Micro Average Recall avg(TPi/(TPi+FNi))= " <<  100*recall << " (%)" << endl;
	cout << "Micro Average F1 avg(2*precision*recall/(precision+recall))= " << 100*F1<< " (%)" << endl;
	cout << "Micro Average Accuracy avg((TPi+TNi)/(TPi+FPi+TNi+FNi))= " <<  100*accuracy << " (%)" << endl;
	cout << "Micro Average ErrorRate avg((FPi+FNi)/TPi+FPi+TNi+FNi))= " <<  100*errorRate << " (%)" << endl;
	cout << "Micro Average Sensitivity = Recall (TPi/(TPi+FNi))=" <<  100*sensitivity << "(%)" << endl;
	cout << "Micro Average Specificity (TNi/(TNi+FPi))=" <<  100*specificity << "(%)" << endl;
	cout << "Micro Average 1 - Specificity (1 - TNi/(TNi+FPi))=" <<  100-100*specificity << "(%)" << endl;
	cout << "AlignmentScore (edition cost): " << alignmentScore << endl;
	cout << "RelativeAlignmentScore (alignment_score/length_GT): " << relativeAlignmentScore << endl;
	cout << "RelativeLatency (Sum_latency/nb_matchs, in #frames): " << relativeLatency << endl;
	cout << "MeanMatchDuration: " << meanDuration << endl;

	cout << endl;
	cout << "Confusion Matrix" << endl;
	for(unsigned int i=beg; i<nbLabs; i++){
		for(unsigned int j=beg; j<nbLabs; j++){
			cout << confusionMatrix[i][j] << " ";
		}
		cout << endl;
	}
}

void EvalTSS::evaluateMetricsFromConfusionMatrix_1vsAll(int beg, int ilab){
	TP=0, FP=0, TN=0, FN=0; N=0;
	int Ni=0;

		TP=confusionMatrix[ilab][ilab];
		for(unsigned int j=beg; j<nbLabs; j++){
			for(unsigned int i=0; i<nbLabs; i++)
	            N+=confusionMatrix[i][j];
			if(ilab!=j)
			   FP+=confusionMatrix[ilab][j];
		}

	nFN=0;
	if(beg==1){

			FN=confusionMatrix[0][ilab];
			nFN=confusionMatrix[0][ilab];
		}

	TN=N-TP-FP-FN;
    	Ni=TP+FN;
	precision=0;
	recall=0;
	accuracy=0;
	specificity=0;
	sensitivity=0;
	errorRate=0;
	globalErrorRate=0;
	F1=0;
	double p,r;


	precision=(float)(TP+1e-10)/(float)(TP+FP+1e-10);
	recall=(float)(TP+1e-10)/(float)(TP+FN+1e-10);
	specificity=(float)TN/(float)(TN+FP+1e-10);
	accuracy=(float)(TP+TN)/(float)(TP+FP+FN+TN+1e-10);
	errorRate+=(float)(FP+FN)/(float)(TP+FP+FN+TN+1e-10);
	F1=2*precision*recall/(precision+recall+1e-10);
	globalErrorRate=TP;

	sensitivity=recall;

	globalErrorRate=errorRate;

}

void EvalTSS::evaluateMetricsFromConfusionMatrix(){
	evaluateMetricsFromConfusionMatrix(0);
}
void EvalTSS::evaluateMetricsFromConfusionMatrix(int beg){
	TP=0, FP=0, TN=0, FN=0; N=0;
    int tp[nbLabs], tn[nbLabs], fp[nbLabs], fn[nbLabs];
cout << beg << " " << nbLabs << endl;
	for(unsigned int i=beg; i<nbLabs; i++){
		tp[i]=confusionMatrix[i][i];
		fp[i]=0; tn[i]=0;
		for(unsigned int j=beg; j<nbLabs; j++){
	        N+=confusionMatrix[i][j];
			if(i!=j)
			   fp[i]+=confusionMatrix[i][j];
		}
	}

	for(unsigned int j=beg; j<nbLabs; j++){
		fn[j]=0;
		for(unsigned int i=beg; i<nbLabs; i++){
			if(i!=j)
			   fn[j]+=confusionMatrix[i][j];
		}
		tn[j]=N-tp[j]-fp[j]-fn[j];
	}
	nFN=0;
	if(beg==1){
		for(unsigned int j=beg; j<nbLabs; j++){
			fn[j]+=confusionMatrix[0][j];
			fp[j]+=confusionMatrix[j][0];
			nFN+=confusionMatrix[0][j];
		}

	}

	precision=0;
	recall=0;
	accuracy=0;
	specificity=0;
	sensitivity=0;
	errorRate=0;
	globalErrorRate=0;
	F1=0;
	double p,r;
	for(unsigned int i=beg; i<nbLabs; i++){
		TP+=tp[i]; TN+=tn[i]; FP+=fp[i]; FN+=fn[i];

		p=(float)(tp[i]+1e-10)/(float)(tp[i]+fp[i]+1e-10);
		precision+=p;
		r=(float)(tp[i]+1e-10)/(float)(tp[i]+fn[i]+1e-10);
		recall+=r;
		specificity+=(float)tn[i]/(float)(tn[i]+fp[i]+1e-10);
		accuracy+=(float)(tp[i]+tn[i])/(float)(tp[i]+fp[i]+fn[i]+tn[i]+1e-10);
		errorRate+=(float)(fp[i]+fn[i])/(float)(tp[i]+fp[i]+fn[i]+tn[i]+1e-10);
		F1+=2*p*r/(p+r+1e-10);
		globalErrorRate+=(float)(tp[i]);
	}
	double nbl=nbLabs-beg;
	precision/=nbl;
	recall/=nbl;
	accuracy/=nbl;
	errorRate/=nbl;
	F1=2*precision*recall/(precision+recall+1e-10);;
	sensitivity=recall;
	specificity=specificity/nbl;
	globalErrorRate=(float)((N+nFN)-globalErrorRate)/(float)(N+nFN);

}

void EvalTSS::initializeConfusionMatrix(int N){
	nbLabs=N;
	confusionMatrix.resize(N, vector<int>(N, 0.0));
}

void EvalTSS::saveFileConfMat(std::string fname){
	std::ofstream fout;
	fout.open(fname.c_str(), ios::out);
	fout << confusionMatrix.size() << endl;
	for(unsigned int i=0; i<confusionMatrix.size(); i++){
		for(unsigned int j=0; j< confusionMatrix[i].size(); j++){
			fout << confusionMatrix[i][j] << " ";
		}
		fout << endl;
	}
	fout.close();
}

void EvalTSS::readFileConfMat(std::string fname){
	std::ifstream fin;
	int val;
	fin.open(fname.c_str(), ios::in);
	fin >> nbLabs;
	initializeConfusionMatrix(nbLabs);
	for(unsigned int i=0; i<nbLabs; i++){
		for(unsigned int j=0; j< nbLabs; j++){
			fin >>val;
			confusionMatrix[i][j]=val;
		}
	}
	fin.close();
}



void EvalTSS::saveFileA(std::string fname, bool append, int nTEST){
	std::ofstream fout;
	bool entete=false;
	if(!std::ifstream(fname) || !append){
		entete=true;
	}
	fout.close();
	if(append)
	   fout.open(fname.c_str(), ios::app);
	else
		fout.open(fname.c_str(), ios::out);

	if(entete)
		fout << "nTEST MERR MACC maPREC maREC maF1 maACC maERR maSENS maSPEC #RPT relRPT AlScr RAlScr relLatcy meanDur" << endl;

	fout << nTEST << " ";
	fout << globalErrorRate << " ";
	fout << 1-globalErrorRate << " ";
	fout << precision << " ";
	fout << recall << " ";
	fout << F1 << " ";
	fout << accuracy << " ";
	fout << errorRate << " ";
	fout << sensitivity << " ";
	fout << specificity << " ";
	fout << RPT << " ";
	fout << relativeRPT << " ";
	fout << alignmentScore << " ";
	fout << relativeAlignmentScore << " ";
	fout << relativeLatency << " ";
	fout << meanDuration << endl;

    fout.close();
}


void EvalTSS::saveFile(std::string fname, bool append, int nTEST){
	std::ofstream fout;
	bool entete=false;
	if(!std::ifstream(fname) || !append){
		entete=true;
	}
	fout.close();
	if(append)
	   fout.open(fname.c_str(), ios::app);
	else
		fout.open(fname.c_str(), ios::out);

	if(entete)
		fout << "nTEST L-GT TP FP TN FN RPT PREC REC F1 ACC SENS SPEC" << endl;

	precision=(float)TP/(float)(TP+FP);
	recall=(float)TP/(float)(TP+FN);
	accuracy=(float)(TP+TN)/(TP+FP+TN+FN);
	F1=2.0*recall*precision/(double)(recall+precision);

	specificity=(float)TN/(float)(TN+FP+1e-10);
	sensitivity=(float)TP/(float)(TP+FN);

	fout << nTEST << " " << Lref << " " << TP << " " <<FP << " " <<TN<< " " << FN<< " " << RPT<< " " << precision<< " " << recall<< " " << F1<< " " << accuracy<< " " << sensitivity<< " " << specificity << endl;
    fout.close();
}

std::vector<int> EvalTSS::loadLabels(std::string &filename, int skip){
	std::vector<int> res={};
	   std::string line;
	   int label;
	   std::ifstream infile (filename.c_str(), std::ifstream::in);
	   if (!infile.good()) {
	        std::cerr << "EvalTSS.cpp::loadLabels, file not found: " << filename << endl;
	        exit(0);
	   }
	   int i=0;
	   while (std::getline(infile, line)) {
		   if(i>=skip){
	         auto &&sline = std::stringstream{line};
	         if(sline >> label){
	        	 res.push_back(label);
	         }
		   }
		   i++;
	   }
	   return res;
}

void EvalTSS::saveStaticLabels(std::string fname, std::vector<int> labs){
	   std::ofstream fout;
	   fout.open(fname.c_str(), ios::out);
	   if (!fout.good()) {
	        std::cerr << "EvalTSS.cpp::saveStaticLabels, could not open file: " << fname << endl;
	        exit(0);
	   }
	   for(unsigned int i=0; i< labs.size(); i++){
		   fout << labs[i]+1 << endl;
	   }
	   fout.close();
}


void EvalTSS::evalDP(std::string &fileRef, std::string &filePred){
	TSS sref=TSS::loadFile(fileRef);
	TSS spred=TSS::loadFile(filePred);
	SeqAlignWithBacktrace(sref, spred);
	evaluateMetricsFromConfusionMatrix(1);
}

void EvalTSS::evalIsolatedPatterns(std::string &fileRef, std::string &filePred){
	std::vector<int> labRef = loadLabels(fileRef,0);
	std::vector<int> labPred = loadLabels(filePred,0);
	int nbLabs=countDistinctSymbols(labRef);
	cout << ">>>>>>>>>>>>>>>> NLABS=" << nbLabs << " " << labRef.size() << endl;
	initializeConfusionMatrix(nbLabs+1);
	for (unsigned int i=0; i< labRef.size(); i++){
		confusionMatrix[labPred[i]][labRef[i]]++;
	}
	evaluateMetricsFromConfusionMatrix();
	std::string fname="static_predictedLabs.txt";
	saveStaticLabels(fname, labPred);
}

EvalTSS::EvalTSS(){
	TP=0; TN=0; FP=0; FN=0; RPT=0; nbLabs=0; Lref=0;
	sensitivity=0; specificity=0; accuracy=0; recall=0; precision=0;
	alignmentScore=0; relativeLatency=0; meanDuration=0;
}
EvalTSS::EvalTSS(int N){
	TP=0; TN=0; FP=0; FN=0; RPT=0; nbLabs=0; Lref=0;
	sensitivity=0; specificity=0; accuracy=0; recall=0; precision=0;
	alignmentScore=0; relativeLatency=0; meanDuration=0;
	confusionMatrix.resize(N, vector<int>(N));
}


