all:
	g++ -g -c Symbol.cpp TSS.cpp EvalTSS.cpp -std=c++11 
	g++ -g evTSS.cpp -o evTSS Symbol.o TSS.o EvalTSS.o -std=c++11

clean:	
	rm *.o evTSS 
