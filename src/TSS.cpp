/*************************************************************************
TSS.cpp (Time Stamped TSS)
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
#include "TSS.h"

using namespace std;

TSS::TSS(){
	nbCat=0;
	data={};
}

TSS::TSS(const std::vector<Symbol> &vs){
	data=vs;
}

std::map<int,int> TSS::getLabels(){
    unsigned int LL=data.size();
    std::map<int,int> tab;
    int l=1;
    for (unsigned int i=0; i<LL; i++){
        if(tab.find(data[i].lab) == tab.end()){
            tab[(int)data[i].lab]=l;
            l++;
        }
    }
    return tab;
}




void TSS::saveFile(std::string file){
ofstream fout;
fout.open(file);
for(unsigned int j=0; j < data.size(); j++){
     fout << data[j].lab << " " << data[j].beg << " " << data[j].end << endl;
 }
fout.close();
}

void TSS::print(){
for(unsigned int j=0; j < data.size(); j++){
     cout << data[j].lab << " " << data[j].beg << " " << data[j].end << endl;
 }
}

TSS TSS::loadFile(std::string &filename){
    auto result = std::vector<Symbol>{};

    std::string line;
    std::ifstream infile (filename.c_str(), std::ifstream::in);
    if (!infile.good()) {
        std::cerr << "file not found: " << filename;
    }
    while (std::getline(infile, line)) {
        auto &&sline = std::stringstream{line};
        Symbol s;
        while (sline >> s.lab >> s.beg >> s.end) {
            result.push_back(s);
        }
    }

    result.shrink_to_fit();
    infile.close();
    return TSS(result);
}

TSS TSS::loadFileSVMPrediction(std::string &filename){
    TSS result;
    result.nbCat=0;
    std::string line, st;
    std::ifstream infile (filename.c_str(), std::ifstream::in);
    if (!infile.good()) {
        std::cerr << "TSS.cpp:75, file not found: " << filename << endl;
        exit(0);
    }
    std::getline(infile, line);
    auto &&enteleLine = std::stringstream{line};
    enteleLine >> st;
    double d;
    while (enteleLine >> d) {
        result.nbCat++;
    }
    int n=0;
    while (std::getline(infile, line)) {
        auto &&sline = std::stringstream{line};
        Symbol s;
        sline >> s.lab; s.beg=n; s.end=n;
        while (sline >> d) {
            s.prob.push_back(d);
        }
        s.prob.shrink_to_fit();
        result.data.push_back(s);
    }

    result.data.shrink_to_fit();
    infile.close();
    return TSS(result);
}


