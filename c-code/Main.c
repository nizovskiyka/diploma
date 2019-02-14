#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>
#include <float.h>


#define solveStep 0.1
#define solveLength 1 //solving on (0,1) interval
#define numSteps 10


/* run this program using the console pauser or add your own getch, system("pause") or input loop */

struct Neuron {
	double moc;
	double mcb;
	double mba;
	double Joc;
	double Jcb;
	double Jba;
	double weight;
	double** eq_sol;
};

//correct
struct Neuron initNeuron(double moc, double mcb, double mba, double Joc, double Jcb, double Jba){
	double sol_res[numSteps][6];
	solveEquations(moc,mcb,mba,Joc,Jcb,Jba,sol_res,solveStep, numSteps);
	struct Neuron neuron = { moc,mcb,mba,Joc,Jcb,Jba,1,sol_res}; //ok?
	return neuron;	
}

//check if correct
struct Neuron* createList(double approx, double interv, double moc_exp, double mcb_exp, double mba_exp, double Joc_exp, double Jcb_exp, double Jba_exp){
	int listLen = 0;//explicit calculation of list length
	while(interv>approx){
		interv-=approx;
		listLen++;
	}
	listLen++;
	listLen = pow(listLen, 6);
	printf("%d\n", listLen);
	
	struct Neuron* neurList;
	neurList = (struct Neuron*)malloc(sizeof(struct Neuron)*listLen);
	
	double moc_max = moc_exp+interv/2;
	double mcb_max = mcb_exp+interv/2;
	double mba_max = mba_exp+interv/2;
	double Joc_max = Joc_exp+interv/2;
	double Jcb_max = Jcb_exp+interv/2;
	double Jba_max = Jba_exp+interv/2;
	
	double moc_st = moc_exp-interv/2;
	double mcb_st = mcb_exp-interv/2;
	double mba_st = mba_exp-interv/2;
	double Joc_st = Joc_exp-interv/2;
	double Jcb_st = Jcb_exp-interv/2;
	double Jba_st = Jba_exp-interv/2;
	
	double moc_min = moc_exp-interv/2;
	double mcb_min = mcb_exp-interv/2;
	double mba_min = mba_exp-interv/2;
	double Joc_min = Joc_exp-interv/2;
	double Jcb_min = Jcb_exp-interv/2;
	double Jba_min = Jba_exp-interv/2;
	
	int i=0;
	
	while(moc_st<moc_max){
		mcb_st=mcb_min;
		while(mcb_st<mcb_max){
			mba_st=mba_min;
			while(mba_st<mba_max){
				Joc_st = Joc_min;
				while(Joc_st<Joc_max){
					Jcb_st=Jcb_min;
					while(Jcb_st<Jcb_max){
						Jba_st=Jba_min;
						while(Jba_st<Jba_max){
							neurList[i]=initNeuron(moc_st,mcb_st,mba_st,Joc_st,Jcb_st,Jba_st);
							i++;
							Jba_st+=approx;
						}
						Jcb_st+=approx;
					}
					Joc_st+=approx;
				}
				mba_st+=approx;
			}
			mcb_st+=approx;
		}
		moc_st+=approx;
	}
	
	return neurList;
}



//check if correct
double transfer(double F){
	if(abs(F)<2*DBL_EPSILON){
		return 1;
	}
	else{
		return ((2*atan(1/F))/M_PI);
	}
}

//code after diff eq method realisation
double checkFunc(struct Neuron neur, struct Neuron neur_rl){
	double neur_res[numSteps][6] = neur.eq_sol;
	double neur_rl_res[numSteps][6] = neur_rl.eq_sol;
	double diff_res[numSteps][6]   //neuron diff array
	int i,j;
	double max[6] = {0,0,0,0,0,0};
	for(i=0;i<6;i++){
		for(j=0;j<numSteps;j++){
			diff_res[j][i] = fabs(neur_rl_res[j][i]-neur_res[j][i])   //difference between neur and neur_rl
			if(diff_res[j][i]>max[i]){    //here searching max among array
				max[i]=diff_res[j][i];
			}
		}
	}
	return pow(max[3],2) + pow(max[4],2) + pow(max[5],2); //return sum of y'
}

//check if correct
struct Neuron getMax(struct Neuron* neurList){
	struct Neuron maxNeuron = neurList[0];
	int length = sizeof(neurList)/sizeof(neurList[0]);//right?
	int i;
	for(i=0;i<length;i++){
		if(neurList[i].weight>maxNeuron.weight){
			maxNeuron = neurList[i];
		}
	}
	return maxNeuron;
}

struct Neuron network(struct Neuron* neurList, struct Neuron neuronRl){
	int i;
	int length = sizeof(neurList)/sizeof(neurList[0]); //right?
	struct Neuron max; //storage size of isn't known
	double check;
	double weightCoeff;
	for(i=0;i<length;i++){
		check = checkFunc(neurList[i],neuronRl);
		weightCoeff = transfer(check);
		neurList[i].weight*=weightCoeff;
	}
	max = getMax(neurList);
	return max;
}

int main(int argc, char *argv[]) {
	
	double approx = 0.2;
	double interv = 0.7;

	
	double human_weight=85;
	
	double moc_exp = 0.1221*human_weight;
    double mcb_exp = 0.0465*human_weight;
    double mba_exp = 0.0146*human_weight;
    double Joc_exp = 0.1137+moc_exp*(0.44*0.44);
    double Jcb_exp = 0.0391+mcb_exp*(0.42*0.42);
    double Jba_exp = 0.0034+mba_exp*(0.44*0.44);
    //real values (to be deleted)
    double moc_rl = 10.37;
    double mcb_rl = 3.9;
    double mba_rl = 1.241;
    double Joc_rl = 2.12;
    double Jcb_rl = 0.736;
    double Jba_rl = 0.24;
    
    struct Neuron* neurList;
    struct Neuron neurRl;
    struct Neuron res;
    
    double solResult[6][numSteps];//maybe won't work
    
    
    neurList = createList(approx,interv,moc_exp,mcb_exp,mba_exp,Joc_exp,Jcb_exp,Jba_exp);
    neurRl = initNeuron(moc_rl,mcb_rl,mba_rl,Joc_rl,Jcb_rl,Jba_rl);
    res = network(neurList,neurRl);//correct?
    
    //RESULT
    
//	printf("result\n");
//    printf("moc: %f\n", moc);
//    printf("mcb: %f\n", mcb);
//    printf("mba: %f\n", mba);
//    printf("Joc: %f\n", Joc);
//    printf("Jcb: %f\n", Jcb);
//    printf("Jba: %f\n", Jba);
    
    
	return 0;
}

