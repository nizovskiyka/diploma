#include <math.h>
#include <stdio.h>

//ode system
double ddalp(double dalp){
	return dalp;
}
double ddbet(double dbet){
	return dbet;
}
double ddgam(double dgam){
	return dgam;
}
double dalp(){
	return 0;
}
double dbet(){
	return 0;
}
double dgam(){
	return 0;
}

double* solveEquations(double moc, double mcb, double Joc, double Jcb, double Jba, double step, double length){

	int i;
	
	// here write initial conditions	
	double dalp_0=0;
	double dbet_0=0;
	double dgam_0=0;
	double alp_0=0;
	double bet_0=0;
	double gam_0=0;
	
	double y1,y2,y3,y4,y5,y6; //function
	double k11,k12,k13,k14,k15,k16; //k1 vector
	double k21,k22,k23,k24,k25,k26; //k2 vector
	double k31,k32,k33,k34,k35,k36; //k3 vector
	double k41,k42,k43,k44,k45,k46; //k4 vector
	
	//init cond set
	y1=dalp_0;
	y2=dbet_0;
	y3=dgam_0;
	y4=alp_0;
	y5=bet_0;
	y6=gam_0;
	
	int numSteps=0;
	while(length<step){
		length-=step;
		numSteps++;
	}

	//now numSteps is the length of the array
	//must be multi-dimentional array
	double result[6][numSteps];
	
	
	//here calculation algorithm
	
	for(i=0;i<numSteps;i++){
		
		//calculating k-vectors
		k11=ddalp(y1);
		k12=ddbet(y2);
		k13=ddgam(y3);
		k14=dalp();
		k15=dbet();
		k16=dgam();
		
		k21=ddalp(y1+step*k11/2);
		k22=ddbet(y2+step*k12/2);
		k23=ddgam(y3+step*k13/2);
		k24=dalp();
		k25=dbet();
		k26=dgam();
		
		k31=ddalp(y1+step*k21/2);
		k32=ddbet(y2+step*k22/2);
		k33=ddgam(y3+step*k23/2);
		k34=dalp();
		k35=dbet();
		k36=dgam();
		
		k41=ddalp(y1+step*k31);
		k42=ddbet(y2+step*k32);
		k43=ddgam(y3+step*k33);
		k44=dalp();
		k45=dbet();
		k46=dgam();
		
		//saving y
		/*
		result_y1[i]=y1;
		result_y2[i]=y2;
		result_y3[i]=y3;
		result_y4[i]=y4;
		result_y5[i]=y5;
		result_y6[i]=y6;
		*/
		result[1][i] = y1;
		result[2][i] = y2;
		result[3][i] = y3;
		result[4][i] = y4;
		result[5][i] = y5;
		result[6][i] = y6;
		
		//calculating next y
		y1 = y1+step*(k11+k21+k31+k41)/6;
		y2 = y2+step*(k12+k22+k32+k42)/6;
		y3 = y3+step*(k13+k23+k33+k43)/6;
		y4 = y4+step*(k14+k24+k34+k44)/6;
		y5 = y5+step*(k15+k25+k35+k45)/6;
		y6 = y6+step*(k16+k26+k36+k46)/6;
	}
	
	return result;
}
