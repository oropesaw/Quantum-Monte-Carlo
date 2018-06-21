// Variational Monte Carlo for the harmonic oscillator

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <ctime>
#include "random.c"
#include "ran2.c"

using namespace std;

int N;                         // number of walkers
double *x;                     // walker positions
double delta;                  // step size
double eSum;                   // accumulator to find energy
double eSqdSum;                // accumulator to find fluctuations in E
double xMin = -10;             // minimum x for histogramming psi^2(x)
double xMax = +10;             // maximum x for histogramming psi^2(x)
double dx = 0.1;               // psi^2(x) histogram step size
double *psiSqd;                // psi^2(x) histogram
int nPsiSqd;                   // size of array
long seed = -time(NULL);

void zeroAccumulators() {
    eSum = eSqdSum = 0;
    for (int i = 0; i < nPsiSqd; i++)
        psiSqd[i] = 0;
}
void initialize() {

    x = new double [N];
    for (int i = 0; i < N; i++)
        x[i] = ran2(&seed) - 0.5;
    delta = 1;

    nPsiSqd = int((xMax - xMin) / dx);
    psiSqd = new double [nPsiSqd];

    zeroAccumulators();
}

double alpha;                  // trial function is exp(-alpha*x^2)

double p(double xTrial, double x) {

    // compute the ratio of rho(xTrial) / rho(x)
    return exp(- 2 * alpha * (xTrial*xTrial - x*x));
}

double eLocal(double x) {

    // compute the local energy
    return alpha + x * x * (0.5 - 2 * alpha * alpha);
}

int nAccept;                   // accumulator for number of accepted steps

void MetropolisStep() {

    // choose a walker at random
    int n = int(ran2(&seed) * N);

    // make a trial move
    double xTrial = x[n] + delta * gasdev(&seed);

    // Metropolis test
    if (p(xTrial, x[n]) > ran2(&seed)) {
        x[n] = xTrial;
        ++nAccept;
    }

    // accumulate energy and wave function
    double e = eLocal(x[n]);
    eSum += e;
    eSqdSum += e * e;
    int i = int((x[n] - xMin) / dx);
    if (i >= 0 && i < nPsiSqd)
        psiSqd[i] += 1;
}

void oneMonteCarloStep() {

    // perform N Metropolis steps
    for (int i = 0; i < N; i++) {
        MetropolisStep();
    }
}

int main() {

   
    N = 10000;
    
    alpha = 0.2;
    
    int MCSteps = 300;
    

    ofstream file("psiSqd.data");
    while(alpha < 1.0){
        initialize();

        // perform 20% of MCSteps as thermalization steps
        // and adjust step size so acceptance ratio ~50%
        int thermSteps = int(0.2 * MCSteps);
        int adjustInterval = int(0.1 * thermSteps) + 1;
        nAccept = 0;
        //cout << " Performing " << thermSteps << " thermalization steps ..."
        //   << flush;
        for (int i = 0; i < thermSteps; i++) {
            oneMonteCarloStep();
            if ((i+1) % adjustInterval == 0) {
                delta *= nAccept / (0.5 * N * adjustInterval);
                nAccept = 0;
            }
        }
        //cout << "\n Adjusted Gaussian step size = " << delta << endl;

        // production steps
        zeroAccumulators();
        nAccept = 0;
        // cout << " Performing " << MCSteps << " production steps ..." << flush;
        for (int i = 0; i < MCSteps; i++)
            oneMonteCarloStep();

    
        double eAve = eSum / double(N) / MCSteps;
        double eVar = eSqdSum / double(N) / MCSteps - eAve * eAve;
        double error = sqrt(eVar) / sqrt(double(N) * MCSteps);
        file << eAve << "\t" << error << "\t" << eVar << "\t" << alpha << endl;
    
        alpha += 0.01;
    }

    file.close();
    
}