#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <ctime>
#include "random.c"
#include "ran2.c"

using namespace std;


const int DIM = 6;             // Dimensión del espacio
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
double rMax = 6;              // máx valor para la función de onda
const int NPSI = 100;          // número de bins del histograma que crea la función de onda.
double psi[NPSI][NPSI];        // histograma de función de ondas



double V(double *r){
    double r_1_mod = 0, r_2_mod = 0, r_12_mod = 0;

    for(int d = 0; d < 3; d++){
        r_1_mod += r[d] * r[d];
        r_2_mod += r[d + 3] * r[d + 3];
        r_12_mod += (r[d] - r[d + 3])*(r[d] - r[d + 3]);
    }

    r_1_mod = sqrt(r_1_mod);
    r_2_mod = sqrt(r_2_mod);
    r_12_mod = sqrt(r_12_mod);

    double pot = 1.0 / r_12_mod - 2.0 / r_2_mod - 2.0 / r_1_mod;
    printf("%lf \n", pot);
    return pot;
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



void zeroAccumulators() {
    ESum = ESqdSum = 0;
    for (int i = 0; i < NPSI; i++)
        for (int j = 0; j < NPSI; j++)
            psi[i][j] = 0;
}


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

   
    for (int d = 0; d < DIM; d++)
        r[n][d] += gasdev(&seed) * sqrt(dt);

    
    double q = exp(- dt * (V(r[n]) - E_T));
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

   
    int N_0 = N;
    for (int n = 0; n < N_0; n++)
        oneMonteCarloStep(n);

   
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

    
    E_T += log(N_T / double(N)) / 10;

   
    ESum += E_T;
    ESqdSum += E_T * E_T;
    
    for (int n = 0; n < N; n++) {
        double rSqd_1 = 0;
        double rSqd_2 = 0;
        for (int d = 0; d < DIM/2; d++){
            rSqd_1 = r[n][d] * r[n][d];
            rSqd_2 = r[n][d + DIM/2] * r[n][d + DIM/2];
        }
        int i = int(sqrt(rSqd_1) / rMax* NPSI);
        int j = int(sqrt(rSqd_2) / rMax* NPSI);
        
        if (i < NPSI && j < NPSI)
            psi[i][j] += 1;
    }
}

int main() {

    N_T = 1000;
    cout << " Enter time step dt: ";
    cin >> dt;
    int timeSteps = 4000;

    initialize();

   
    int thermSteps = int(0.2 * timeSteps);
    for (int i = 0; i < thermSteps; i++)
        oneTimeStep();

    //ofstream folder("time.txt");
    zeroAccumulators();
    for (int i = 0; i < timeSteps; i++) {
        oneTimeStep();
       // folder  << i* dt<< "\t" << N << "\t" << E_T << endl;
    }
    //folder.close();
   
    double EAve = ESum / timeSteps;
    double EVar = ESqdSum / timeSteps - EAve * EAve;
    
    cout << " <E> = " << EAve << " +/- " << sqrt(EVar / timeSteps) << endl;
    cout << " <E^2> - <E>^2 = " << EVar << endl;
   
    double psiNorm = 0, psiExactNorm = 0;
    double dr = rMax / NPSI;
    
    for (int i = 0; i < NPSI; i++)
        for (int j = 0; j < NPSI; j++){
        psiNorm += dr * dr * psi[i][j] * psi[i][j];
    }

    psiNorm = sqrt(psiNorm);
   
    ofstream file("psi.data");
    
    for (int i = 0; i < NPSI; i++) 
        for (int j = 0; j < NPSI; j++){
        double r_1 = i * dr;
        double r_2 = j * dr;
        file << r_1 << "\t" << r_2 << '\t' << psi[i][j] / psiNorm << endl;
    }
    file.close();
}
