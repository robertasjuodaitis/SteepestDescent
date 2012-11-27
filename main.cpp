//
//  main.cpp
//  SteepestDescent
//
//  Created by Remigijus Paulavicius on 11/20/12.
//  Copyright (c) 2012 Remigijus Paulavicius. All rights reserved.
//

#include <iostream>
#include <cmath>
#include <cstdlib> // rand(), srand()
#define N 2
using namespace std;

// Vektoriaus begalines (max) normos funkcijos deklaracija 
double Vector_Max_Norm(double v[], int n);

// Greiciausio nusileidimo (angl. Steepest Descent) metodo deklaracija
int  Steepest_Descent(double (*f)(double *), void (*df)(double *, double *),
     int (*stopping_rule)(double*, double, double*, double, double*, int, int), 
                          double a[], double *fa, double *dfa, double cutoff,
						double cutoff_scale_factor, double tolerance, int n);

// Generuoja atsitiktini realu skaiciu tarp dLow and dHigh
double GetRandomNumber(double dLow, double dHigh){
    return static_cast<double>(rand())/RAND_MAX*(dHigh-dLow) + dLow;
}

// Apskaiciuoja Six-hump Camel Back funkcijos reiksme taske x
double SixHumpCamelBack(double *x){
    return (4-2.1*x[0]*x[0]+x[0]*x[0]*x[0]*x[0]/3)*x[0]*x[0] + x[0]*x[1] +
    (-4+4*x[1]*x[1])*x[1]*x[1];
}
// Apskaiciuoja Six-hump Camel Back gradiento reiksme taske x
void SixHumpCamelBackGradient(double *x, double *fGrad){
    fGrad[0] = 8*x[0]-8.4*x[0]*x[0]*x[0]+2*x[0]*x[0]*x[0]*x[0]*x[0]+x[1];
    fGrad[1] = x[0]-8*x[1]+16*x[1]*x[1]*x[1];
}

// Algoritmo sustojimo salyga kontroliuojanti funkcija
int StoppingRule(double* a, double fa, double* x, double fx, double* dfa, int
iteration, int n){
	double fEps = abs(fx - fa); // Funkcijos reiksmiu skirtumas
	double xa[n];
	for(int i = 0; i < n; ++i) xa[i] = x[i]-a[i];
	double xEps = Vector_Max_Norm(xa, 2); // Argumento skirtumo norma
	double dfaEps = Vector_Max_Norm(dfa, 2); // Gradiento norma
	if(iteration > 3) 
		return -6;
	else
		return 0;
}

void RandomSearch(double *a);

int main(int argc, const char * argv[])
{
	double region[] = {-1.9, 1.9, -1.1, 1.1};
	
	//tik randomui
    double a[6];
	RandomSearch(a);

//kol kas isjungiu
//    double a[N] = {0.0, 1.0}; // N-matis Vektorius
	
    
    /*srand(time(0)); // Naudoja vis kita seed'a
	double a[N]; // N-matis Vektorius
	for(int i = 0; i < N; ++i){
        a[i] = GetRandomNumber(region[2*i], region[2*i+1]);
    }*/
	double fa = SixHumpCamelBack(a); // Funkcijos reiksme pradiniame taske a
	double dfa[N]; 
	SixHumpCamelBackGradient(a, dfa); // Funkcijos gradiento reiksme taske a
	double cutoff = 1.0, cutoff_scale_factor = 1.0; // Pap. parametrai
	double tolerance = 0.01;
	int err = Steepest_Descent( SixHumpCamelBack, SixHumpCamelBackGradient, StoppingRule,
    a, &fa, dfa, cutoff, cutoff_scale_factor, tolerance, N);
	switch (err)
	{
		case 0:
			cout << "Success" << endl;
			break;
		case -1:
			cout << "In the line search three points are collinear." << endl; 
			break;
		case -2:
			cout << "In the line search the extremum of the parabola through the three points is a maximum." << endl;
			break;
		case -3:
			cout << "Int the line search the initial points failed to satisfy the condition that x1 < x2 < x3 and fx1 > fx2 < fx3." << endl;
			break;
		case -4:
			cout << "Not enough HEAP memory." << endl;
			break;
		case -5:
			cout << "The gradient evaluated at the initial point vanishes." << endl;
		case -6:
			cout << "Exceed maximal number of iterations." << endl;
		break;
	}
	cout << "Greiciausio nusileidimo (angl. Steepest Descent) metodu" << endl;
	cout << "surastas sprendinys yra:" << endl;
	cout << "xMin = (" << a[0] << ", " << a[1] << ")" << endl;
	cout << "f(xMin) = " << fa << endl;
	system("pause");
	return 0;
	
}

