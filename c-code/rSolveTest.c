#include "rkSolve.h"

int main(int argc, char *argv[]){
	double approx = 0.2;
	double interv = 0.7;
	double solveStep = 0.1;
	double solveLength = 1; //solving on (0,1) interval
	
	double moc_rl = 10.37;
    double mcb_rl = 3.9;
    double mba_rl = 1.241;
    double Joc_rl = 2.12;
    double Jcb_rl = 0.736;
    double Jba_rl = 0.24;
    
    int i,j;
    
    int numSteps=0;
	while(solveLength>solveStep){
		solveLength-=solveStep;
		numSteps++;
	}
	
	double solResult[numSteps][6];//maybe won't work
	
	solveEquations(moc_rl, mcb_rl, mba_rl, Joc_rl, Jcb_rl, Jba_rl, solResult, solveStep, numSteps);
	/*
	for(i=0;i<6;i++){
		for(j=0;j<numSteps;j++){
			printf("%f, ", solResult[j][i]);
		}
		printf("a\n");
	}*/
	return 0;
}
