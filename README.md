
# evTSS Program 

**Pairwise Alignment of Time Stamped Sequences for assessing automatic labelling of stream data given a ground truth**

**evTSS** takes in input two file names (two Time Stamped Sequences (TSS)), the first one containing the timestamped ground truth labels and the second one the timestamped predicted labels. It provides as outputs some alignment assessment metrics, as detailed in <https://hal.archives-ouvertes.fr/hal-01403948v1>

**evTSS** can be usefull to assess pattern detection and recognition in stream as well as automatic labelling tasks, provided some ground truth is available.

**evTSS** is available as C++ code, and has been tested on recent Linux Ubuntu distributions.

## Installation
to compile evTSS, just execute in the src directory
$ make

## Usage
in the src directory run the program as 
$ ./evTSS [GroundTruth_fileName]  [Predicted_fileName] 

**[GroundTruth]** file name of the sequence containing the gound truth symbols
**[Predicted]**  file name of the sequence containing the predicted symbols



E.g.: 

$ ./evTSS ./test/groundtruth.txt ./test/predicted.txt


Each input file contains a sequence of symbols: each row of the input files contains a symbol, each symbol being a triplet *[label, begin\_ index, end\_ index]* where  *label* is in \{0,..,|A|\} and *N* being the size of the Alphabet *A=\{1,2, ... ,|A|\}*, and *begin_index* and *end_index* are integers or floats.

**WARNING:** ** 0;   ;** is an extra symbol used to indicate that **no label** has been emitted (in the ground truth) or detected (in the predicted sequence). 

### Example: Groundtruth.txt

3 0 45
0 46 50
5 51 101
2 102 152
4 153 203
1 204 254

// symbol 1 starts at time index 0 and ends at time index 50
// symbol 3 starts at time index 51 and ends at time index 101
...
// symbol 4 starts at time index 255 and ends at time index 305

### Example: Predicted.txt

0 0 30
0 31 50
0 51 88
5 89 90
0 91 95
5 96 106
2 107 152
0 153 174
2 175 195
0 196 203
1 204 254

// Nothing is detected on time segment [0-20]
// symbol 6 starts at time index 21 and ends at time index 40
// Nothing is detected on time segment [41-42] 
// symbol 6 starts at time index 43 and ends at time index 47
...
// symbol 5 starts at time index 270 and ends at time index 305

## The program outputs the following trace and evaluation metrics 
**GROUND TRUTH	PREDICTION		EVENTS**
5 (1|204-254)			10 (1|204-254)	*correct match type-1*
4 (4|153-203)			9 (0|196-203)	 	*suppress NL in pred type-3 *
4 (4|153-203)			8 (2|175-195)	 	*mismatch (1 FP and 1 FN) type-1*
3 (2|102-152)			7 (0|153-174)	 	*suppress NL in pred type-3 *
3 (2|102-152)			6 (2|107-152)		*correct match type-1*
2 (5|51-101)			5 (5|96-106)	 	*correct match type-1*
1 (0|46-50)			4 (0|91-95)	 	*suppress NL in pred type-3 *
1 (0|46-50)			3 (5|89-90)	 	*suppress in pred (FP) type-3*
1 (0|46-50)			2 (0|51-88)	 	*suppress NL in pred type-3* 
1 (0|46-50)			1 (0|31-50)	 	*NL match type-1*
0 (3|0-45)			0 (0|0-30)	 	*suppress in GT (FN) type-2*
0 (3|0-45)			0 (0|0-30)	 	*suppress NL in pred type-3* 

Macro Error Rate (N-sum_i(FPi+FNi))/N= 40 (%)
Micro Average Precision avg(TPi/(TPi+FPi))= 80 (%)
Micro Average Recall avg(TPi/(TPi+FNi))= 60 (%)
Micro Average F1 avg(2*precision*recall/(precision+recall))= 68.5714 (%)
Micro Average Accuracy avg((TPi+TNi)/(TPi+FPi+TNi+FNi))= 82 (%)
Micro Average ErrorRate avg((FPi+FNi)/TPi+FPi+TNi+FNi))= 18 (%)
Micro Average Sensitivity = Recall (TPi/(TPi+FNi))=60(%)
Micro Average Specificity (TNi/(TNi+FPi))=88.3333(%)
Micro Average 1 - Specificity (1 - TNi/(TNi+FPi))=11.6667(%)
AlignmentScore (edition cost): 4.00909
RelativeAlignmentScore (alignment_score/length_GT): 0.668182
RelativeLatency (Sum_latency/nb_matchs, in #frames): 9.16667
MeanMatchDuration: 35



**Confusion Matrix**
1 0 0 1 0 0 
0 1 0 0 0 0 
0 0 1 0 1 0 
0 0 0 0 0 0 
0 0 0 0 0 0 
1 0 0 0 0 1 




