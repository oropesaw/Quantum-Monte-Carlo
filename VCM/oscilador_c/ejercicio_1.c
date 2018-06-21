#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>
#include <math.h>

#define squaret(A) (A)*(A)
#define N 5000
#define DELTA 4
# define IM1 2147483563
# define IM2 2147483399
# define AM (1.0/IM1)
# define IMM1 (IM1-1)
# define IA1 40014
# define IA2 40692
# define IQ1 53668
# define IQ2 52774
# define IR1 12211
# define IR2 3791
# define NTAB 32
# define NDIV (1+IMM1/NTAB)
# define EPS 1.2e-7
# define RNMX (1.0-EPS)


typedef double(* function)(double, double);

typedef struct system{
    function f, E_coeff;
    long seed;
    double coeff;
    double energy;
}class_system;

class_system *builder_class(function f, function E_coeff, double coeff);
double f(double x, double coeff);
double E_coeff(double x, double coeff);
void metropolis(class_system *sys);
float ran2(long *idum);

void main(){
    class_system *sys;
    double coeff;
    FILE *folder = fopen("energy.txt", "w");
    assert(folder);

    for(coeff = 0.4; coeff <= 0.6; coeff += 0.01){
        sys = builder_class(f, E_coeff, coeff);
        metropolis(sys);
        fprintf(folder,"%lf\t%lf\n", sys -> energy / N, coeff);
        free(sys);
    }

    fclose(folder);
}

class_system *builder_class(function f, function E_coeff, double coeff){
    
    class_system *sys = (class_system *)malloc(sizeof(class_system));
    assert(sys);

    sys -> f = f;
    sys -> E_coeff = E_coeff;
    sys -> coeff = coeff;
    sys -> seed = -time(NULL);
    
    return sys;
}

double f(double x, double coeff){
    return exp(-coeff * squaret(x));
}

double E_coeff(double x, double coeff){
    return coeff + squaret(x) * (0.5 - 2 * squaret(coeff));
}

void metropolis(class_system *sys){
    double x_old = 0;
    double x_new, p_new, p_old, k;
    int i;
    
    sys -> energy = 0.0;
    
    for(i = 0; i < N; i++){
        x_new = x_old + (ran2(&sys -> seed) - 0.5) * DELTA;

        p_new = squaret(sys -> f(x_new, sys -> coeff));
        p_old = squaret(sys -> f(x_old, sys -> coeff));

        k = (double)(p_new / p_old);
        if(k >= ran2(&sys -> seed)){
             sys -> energy += sys -> E_coeff(x_new, sys -> coeff);
             x_old = x_new;
        } else
            sys -> energy += sys -> E_coeff(x_old, sys -> coeff);
    }

}

float ran2(long *idum){
    int j;
    long k;
    static long idum2=123456789;
    static long iy=0;
    static long iv[NTAB];
    float temp;

    if(*idum <= 0){                        //Initialize.
        if(-(*idum) < 1) *idum = 1;              //Be sure to prevent idum = 0.
        else *idum = -(*idum);
            idum2 = (*idum);
        for (j = NTAB + 7; j >= 0; j--){                //Load the shuffle table (after 8 warm-ups).
            k = (*idum) / IQ1;
            *idum = IA1 * (*idum - k * IQ1) - k * IR1;
            if (*idum < 0) *idum += IM1;
                if (j < NTAB) iv[j] = *idum;
        }
        iy=iv[0];
    }

    k = (*idum) / IQ1;                           //Start here when not initializing.
   *idum = IA1 * (*idum - k * IQ1) - k * IR1;           //Compute idum=(IA1*idum) % IM1 without
    if (*idum < 0) *idum += IM1;             //overflows by Schrage’s method.
        k=idum2/IQ2;
    idum2 = IA2 * (idum2 - k * IQ2) - k * IR2;           //Compute idum2=(IA2*idum) % IM2 likewise.
    if (idum2 < 0) idum2 += IM2;
        j = iy / NDIV;                               //Will be in the range 0..NTAB-1.
    iy = iv[j] - idum2;                          //Here idum is shuffled, idum and idum2 are
    iv[j] = *idum;                           //combined to generate output.
    
    if (iy < 1) iy += IMM1;
    if ((temp=AM*iy) > RNMX) return RNMX;    //Because users don’t expect endpoint values.
    else return temp;
}
