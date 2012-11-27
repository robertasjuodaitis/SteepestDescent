#include <iostream>
#include <cstdlib> // rand(), srand()
#include <ctime> // time() - gražina sekundžiu skaiciu nuo 1970 m. sausio 1 d.
using namespace std;
static int n = 2; // Uzdavinio dimensija

// Generuoja atsitiktini realu skaiciu tarp dLow and dHigh
double GetRandomNumber(double dLow, double dHigh){
    return static_cast<double>(rand())/RAND_MAX*(dHigh-dLow) + dLow;
}
// Apskaiciuoja Six-hump Camel Back funkcijos reiksme taske x
double SixHumpCamelBack(double *x){
    return (4-2.1*x[0]*x[0]+x[0]*x[0]*x[0]*x[0]/3)*x[0]*x[0] + x[0]*x[1] +
    (-4+4*x[1]*x[1])*x[1]*x[1];
}
// Isrinkimo paieskos algoritmas
void Selection(double *dF, double *dX, int l, int r){
    for (int i = l; i < r; i++){
        int min = i;
        for (int j = i+1; j <= r; j++)
            if (dF[j] < dF[min]) min = j;
        double tF = dF[i]; dF[i] = dF[min]; dF[min] = tF;
        double tX[n];
        for (int ii = 0; ii < n; ii++){
            tX[ii] = dX[i*n+ii];
        }
        for (int ii = 0; ii < n; ii++){
            dX[i*n+ii] = dX[min*n+ii];
        }
        for (int ii = 0; ii < n; ii++){
            dX[min*n+ii] = tX[ii];
        }
    }
}
void RezultatuIsvedimas(int nCount, double dGlobF, double *dX){
    cout << nCount << " band. : \n";
    cout << "globF = " << dGlobF << endl;
    for(int i = 0; i < n; ++i)
        cout << "x[" << i << "] = " << dX[i] << endl;
    cout << endl;
}
void RandomSearch(double *a){
    int N = 1000; // Iteraciju skaitiklis
    int count = 1;
    int stopingRule = 1; // 0: FOR N zingsniu; 1: WHILE > eps
    double x[n];
    double *xAll = new double[N*n];
    double f;
    double *fAll = new double[N];
    double globF = 1000000;
    double globX[n];
    double globMinimum = -1.031628453;
    double region[] = {-1.9, 1.9, -1.1, 1.1};
    double eps = 0.01;
  
    srand(time(0)); // Naudoja vis kita seed'a
    if (stopingRule == 0){
        for (count=0; count < N; count++)
        {
            for(int i = 0; i < n; ++i){
                x[i] = GetRandomNumber(region[2*i], region[2*i+1]);
                xAll[count*n+i] = x[i];
            }
            f = SixHumpCamelBack(&x[0]);
            fAll[count] = f;
            if (f < globF){
                globF = f;
                for(int i = 0; i < n; ++i) globX[i] = x[i];
                RezultatuIsvedimas(count, globF, &globX[0]);
            }
        }
    } else {
        while (globF - globMinimum > eps)
        {
            for(int i = 0; i < n; ++i){
                x[i] = GetRandomNumber(region[2*i], region[2*i+1]);
                xAll[count*n+i] = x[i];
            }
            f = SixHumpCamelBack(&x[0]);
            fAll[count] = f;
            if (f < globF){
                globF = f;
                for(int i = 0; i < n; ++i) globX[i] = x[i];
                RezultatuIsvedimas(count, globF, &globX[0]);
            }
            count++;
            if (count >= N){
                double *xAllTemp = xAll;
                double *fAllTemp = fAll;
                N = 2*N; // Padvigubiname masyvo dydi
                xAll = new double[N*n];
                fAll = new double[N*n];
                for (int i=0; i<count; i++)
                {
                    xAll[i] = xAllTemp[i];
                    fAll[i] = fAllTemp[i];
                }
                // Isvalome atminty (HEAP'a)
                delete[] xAllTemp;
                delete[] fAllTemp;
            }
        }
        
    }
    // Isvalome atminty (HEAP'a)
    delete[] xAll;
    delete[] fAll;    cout << "RANDOM SEARCH REZULTATAS: \nPo ";
    RezultatuIsvedimas(count-1, globF, &globX[0]);
    cout << "Pries rusiavima:" << endl;
    for (int i=0; i<count; i++) {
        RezultatuIsvedimas(i, fAll[i], &xAll[i*n]);
    }
    Selection(&fAll[0], &xAll[0], 0, count-1);
    cout << "Po rusiavimo:" << endl;
    for (int i=0; i<count; i++) {
        RezultatuIsvedimas(i, fAll[i], &xAll[i*n]);
    }
    cout<<"Isvedam masyva\n";
    for (int i=0;i<6;i++){
        a[i]= xAll[i];
        cout<<a[i]<<endl;//<<" "<<a[i++]<<endl;
        
    }
}
