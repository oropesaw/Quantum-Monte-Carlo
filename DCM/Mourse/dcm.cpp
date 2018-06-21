#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <ctime>
#include "random.c"
#include "ran2.c"

using namespace std;


const int DIM = 1;             // Dimensión del espacio
long seed = -time(NULL);

double dt;                     // paso de tiempo a usar
double E_T;                    // energía objetivo

// caminantes random
int N;                         // número actual de caminantes
int N_T;                       // número deseado de caminantes
double **r;                    // x,y,z posición de caminantes
bool *alive;                   // verifica cuantos caminantes está muertos.

// observables
double ESum;                   // acumulador de energía
double ESqdSum;                // acumulador de varianza
double rMin = -5;
double rMax = 15;              // máx valor para la función de onda
const int NPSI = 100;          // número de bins del histograma que crea la función de onda.
double psi[NPSI];              // histograma de función de ondas



// Oscilador armónico en DIM dimensión
double V_1(double *r) {                     
    double rSqd = 0;
    for (int d = 0; d < DIM; d++)
        rSqd += r[d] * r[d];
    return 0.5 * rSqd;
}


// Morse oscilador en DIM dimensión
double V_2(double *r){
    double rSqd = 0;
    for (int d = 0; d < DIM; d++)
        rSqd += r[d];//*r[d];
    return 0.5 * (exp(-2*rSqd) - 2.0*exp(-rSqd));
}


// redimensiona el objeto double ** que apunta a los caminantes
// cada vez que uno de estos pueda morir.
void ensureCapacity(int index) {

    static int maxN = 0;       

    if (index < maxN)          
        return;                

    int oldMaxN = maxN;        
    if (maxN > 0)
        maxN *= 2;             
    else
        maxN = 1;
    if (index > maxN - 1)      
        maxN = index + 1;      

    
    double **rNew = new double* [maxN];
    bool *newAlive = new bool [maxN];
    
    for(int n = 0; n < maxN; n++) {
        rNew[n] = new double [DIM];
        if (n < oldMaxN) {     
            for (int d = 0; d < DIM; d++)
                rNew[n][d] = r[n][d];
            newAlive[n] = alive[n];
            delete [] r[n];    
        }
    }
    delete [] r;               
    r = rNew;                  
    delete [] alive;
    alive = newAlive;
}


// acumulador de ceros
void zeroAccumulators() {
    ESum = ESqdSum = 0;
    for (int i = 0; i < NPSI; i++)
        psi[i] = 0;
}


// inicializa el sistema.
void initialize() {
    N = N_T;                   
    for (int n = 0; n < N; n++) {
        ensureCapacity(n);
        for (int d = 0; d < DIM; d++)
            r[n][d] = ran2(&seed) - 0.5;
        alive[n] = true;
    }
    
    zeroAccumulators();
    E_T = 0;                  
}


void oneMonteCarloStep(int n) {

    // paso difusivo
    for (int d = 0; d < DIM; d++)
        r[n][d] += gasdev(&seed) * sqrt(dt);


    // paso de ramificación
    double q = exp(- dt * (V_2(r[n]) - E_T));
    int survivors = int(q);
    if (q - survivors > ran2(&seed))
        ++survivors;

    

    for (int i = 0; i < survivors - 1; i++) {
        ensureCapacity(N);
        for (int d = 0; d < DIM; d++)
            r[N][d] = r[n][d];
        alive[N] = true;
        ++N;
    }

    
    if (survivors == 0)
        alive[n] = false;
}



void oneTimeStep() {

    // paso difusivo por cada caminante.
    int N_0 = N;
    for (int n = 0; n < N_0; n++)
        oneMonteCarloStep(n);

    // eliminación de los caminantes muertos.
    int newN = 0;
    for (int n = 0; n < N; n++)
    if (alive[n]) {
        if (n != newN) {
            for (int d = 0; d < DIM; d++)
                r[newN][d] = r[n][d];
            alive[newN] = true;
        }
        ++newN;
    }
    N = newN;

    // ajuste de  E_T
    E_T += log(N_T / double(N)) / 10;

    // registro de energía y función de ondas.
    ESum += E_T;
    ESqdSum += E_T * E_T;
    
    for (int n = 0; n < N; n++) {
        double rSqd = 0;
        for (int d = 0; d < DIM; d++)
            rSqd = r[n][d];
        int i = int(rSqd / (rMax - rMin) * NPSI) + int(abs(rMin) / (rMax - rMin) * NPSI);
        
        if (i < NPSI)
            psi[i] += 1;
    }
}


int main() {

    N_T = 300;
    dt = 0.01;
    int timeSteps = 4000;
    
    ofstream file_0("data.txt");

    while(dt <= 0.2){
        initialize();

        int thermSteps = int(0.2 * timeSteps);
        for(int i = 0; i < thermSteps; i++)
            oneTimeStep();


        zeroAccumulators();
        for(int i = 0; i < timeSteps; i++)
            oneTimeStep();

    
        double EAve = ESum / timeSteps;
        double EVar = ESqdSum / timeSteps - EAve * EAve;
    
        file_0 << EAve << "\t" << sqrt(EVar / timeSteps) << "\t" << EVar << "\t" << dt << endl;
        dt += 0.01;
    }

    double psiNorm = 0, psiExactNorm = 0;
    double dr = (rMax - rMin) / NPSI;
    
    for (int i = 0; i < NPSI; i++) {
        double r = i * dr + rMin;
        
        psiNorm += dr * psi[i] * psi[i];
        psiExactNorm +=  dr * exp(-2*exp(- r) - r);
    }

    psiNorm = sqrt(psiNorm);
    psiExactNorm = sqrt(psiExactNorm);
    
    ofstream file_1("psi.data");
    
    for (int i = 0; i < NPSI; i++) {
        double r = i * dr + rMin;
        file_1 << r << '\t' << psi[i] / psiNorm << '\t'
             <<  exp(-exp(- r) - 0.5* r) / psiExactNorm << '\n';
    }
    file_1.close();
}
