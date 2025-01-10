#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <random>
#include <vector>
using namespace std;

random_device rd;
mt19937_64 mt64(rd());

double uniform(){
    uniform_real_distribution<double> get_rand_uni_real(0.,1.);

    return get_rand_uni_real(mt64);
}

uint64_t get_rand() {
    return mt64();
}

#define L 100 //systme size
#define T 2.5 //temperature
#define tau0 1. // typical time
#define Tg 1.
#define n 2.

gamma_distribution<double> gamma_dis(n,Tg);

int main(int argc, char *argv[]){
    FILE *fp;
    char file[100];
    int i,j,k;
    double tau,E;

    snprintf(file,100,"gamma/tau1.txt");
    fp=fopen(file,"w");
    for(i=0;i<L;++i){
        E=gamma_dis(mt64);
        tau=tauc*exp(E/T);

        fprintf(fp,"%lf\n",tau);
    }

    fprintf(fp,"\n#n=%.1f Tg=%.1f T=%.1f",n,Tg,T);

    fclose(fp);
}
