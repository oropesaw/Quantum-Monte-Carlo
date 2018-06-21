#include <iostream>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <ctime>
#include "ran2.c"

using namespace std;

const int DIM = 3;                    // dimensión del espacio.
const int NEL = 2;                    // número de electrones.
int N = 50;                           // número de caminantes.
double (*R)[NEL][DIM];                // cordenadas en espacio de configuraciones.

double beta;                          // parámetro variacional 1.
double alpha;                         // parámetro variacional 2.

double delta;                         // parámetro de aceptación;
long seed = -time(NULL);              // semilla del generador aleatorio.

int N_acp;
double E_sum;
double E_sqd_sum;

void initialize();
double Psi(double *r_1, double *r_2);
double E_loc(double *r_1, double *r_2);
void metropolis(int walker);
void one_monte_carlo_step();
void zero_accumulators();


int main() {
    int MCSteps = 1000;
    ofstream folder("data.txt");

    for(beta = 1.8; beta < 2.1; beta += 0.01){
        for(alpha = 0.0; alpha < 0.6; alpha += 0.01){
            initialize();

            // realiza el 20% de MCSteps como pasos de termalización
            // y ajusta el tamaño del paso, por lo que la tasa de aceptación es ~50%
            int thermSteps = int(0.2 * MCSteps);
            int adjustInterval = int(0.1 * thermSteps) + 1;
            N_acp = 0;
            
            for (int i = 0; i < thermSteps; i++) {
                one_monte_carlo_step();
                if ((i+1) % adjustInterval == 0) {
                    delta *= N_acp / (0.5 * N * adjustInterval);
                    N_acp = 0;
                }
            }
            

            
            zero_accumulators();
            N_acp = 0;
    
            for (int i = 0; i < MCSteps; i++)
            one_monte_carlo_step();
            
            
            // cálculo de la energía. 
            double eAve = E_sum / double(N) / MCSteps;
            double eVar = E_sqd_sum / double(N) / MCSteps - eAve * eAve;
            double error = sqrt(eVar) / sqrt(double(N) * MCSteps);
            folder << eAve << "\t" << error << "\t"<< eVar << "\t" << beta << "\t" << alpha << endl;

        }
    }
    folder.close();
}


// función que inicializa el sistena.
void initialize(){
    R = new double [N][NEL][DIM];
    for(int n = 0; n < N; n++)
        for(int e = 0; e < NEL; e++)
            for(int d = 0; d < DIM; d++)
                R[n][e][d] = ran2(&seed) - 0.5;
    delta = 1;
}


// función de acumulación de ceros.
void zero_accumulators(){
    E_sum = E_sqd_sum = 0;
}


// función de onda del oscilador.
double Psi(double *r_1, double *r_2){
    // valor de la función de onda de prueba para el caminante n.
    double r_1_mod = 0, r_2_mod = 0, r_12_mod = 0;
    
    for(int d = 0; d < 3; d++){
        r_1_mod += r_1[d] * r_1[d];
        r_2_mod += r_2[d] * r_2[d];
        r_12_mod += (r_1[d] - r_2[d])*(r_1[d] - r_2[d]);
    }
    r_1_mod = sqrt(r_1_mod);
    r_2_mod = sqrt(r_2_mod);
    r_12_mod = sqrt(r_12_mod);
    
    double Psi = - beta * (r_1_mod + r_2_mod) + r_12_mod / (2 * (1 + alpha*r_12_mod));
    return exp(Psi);
}


// energía local del oscilador.
double E_loc(double *r_1, double *r_2){
    double r_1_mod = 0, r_2_mod = 0, r_12_mod = 0;
    
    for(int d = 0; d < 3; d++){
        r_1_mod += r_1[d] * r_1[d];
        r_2_mod += r_2[d] * r_2[d];
        r_12_mod += (r_1[d] - r_2[d])*(r_1[d] - r_2[d]);
    }
    r_1_mod = sqrt(r_1_mod);
    r_2_mod = sqrt(r_2_mod);
    r_12_mod = sqrt(r_12_mod);
    
    double dopProd = 0;
    
    for(int d = 0; d < 3; d++){
        dopProd += (r_1[d]*r_2[d]) / (r_1_mod*r_2_mod);
    }
    
    double denom =  (1 + alpha * r_12_mod);
    double denom2 = 2 * r_12_mod * denom * denom;
    double denom3 = r_12_mod * denom * denom;
    double denom4 = denom * denom * denom;
    double denom5 = 4 * denom * denom4;
    
    double e = - beta * beta + (beta - 2) * (1.0 / r_1_mod + 1.0 / r_2_mod) + 1.0 / r_12_mod +
        (beta * (r_1_mod + r_2_mod) * (1 - dopProd)) / denom2 - 1.0 / denom3 + alpha / denom4 - 1.0 / denom5;
    
    return e;
}


// implementación del algoritmo de metrópolis.
void metropolis(int walker){
    //hacer un movimiento de prueba de cada electrón.
    double r_1[3], r_2[3], rTrial1[3], rTrial2[3];

    for (int d = 0; d < 3; d++) {
        r_1[d] = R[walker][0][d];
        rTrial1[d] = r_1[d] + delta * (2 * ran2(&seed) - 1);
        r_2[d] = R[walker][1][d];
        rTrial2[d] = r_2[d] + delta * (2 * ran2(&seed) - 1);
    }

    for (int d = 0; d < 3; d++) {
        r_1[d] = R[walker][0][d];
        rTrial1[d] = r_1[d] + delta * (2 * ran2(&seed) - 1);
        r_2[d] = R[walker][1][d];
        rTrial2[d] = r_2[d] + delta * (2 * ran2(&seed) - 1);
    }

    // Metrópolis test
    double w = Psi(rTrial1, rTrial2) / Psi(r_1, r_2);
    if (ran2(&seed) < w * w) {
        for (int d = 0; d < 3; d++) {
            R[walker][0][d] = r_1[d] = rTrial1[d];
            R[walker][1][d] = r_2[d] = rTrial2[d];
        }
        ++N_acp;
    }

    // acumulador de energía local
    double e = E_loc(r_1, r_2);
    E_sum += e;
    E_sqd_sum += e * e;
}


void one_monte_carlo_step() {
    for (int n = 0; n < N; n++)
    metropolis(n);
}

