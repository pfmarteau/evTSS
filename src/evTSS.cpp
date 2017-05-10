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
#include <iostream>
#include <cfloat>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "EvalTSS.h"


void usage(){
    cout << "./EvalStream  <GroundTruth> <Predicted>" << endl;
    cout << "<GroundTruth>: file name of the sequence containing the gound truth symbols" << endl;
    cout << "<Predicted>: file name of the sequence containing the predicted symbols" << endl;

 }

int main(const int argc, const char **argv) {

    if(argc<3){
        usage();
        return 1;
    }

    string fileRef=argv[1];
    string filePred=argv[2];
    cout << "Evaluate " << filePred << " against " << fileRef << endl;

    EvalTSS ev;

    ev.evalDP(fileRef, filePred);
    ev.print(0);

    return 0;
}


