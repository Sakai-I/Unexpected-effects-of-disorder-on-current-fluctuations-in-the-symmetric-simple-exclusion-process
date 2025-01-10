#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <random>
#include <vector>
using namespace std;

random_device rd;
mt19937_64 mt64(rd());

double uniform(){
    uniform_real_distribution<double> get_rand_uni_real( 0.0, 1.0 );

    // 乱数を生成
    return get_rand_uni_real(mt64);
}

#define L 100 //system size
#define tau0 1. //typical time
#define PI acos(-1.)
#define sigma2 4. 
#define T 1. //temperature

int main(){
    FILE *fp;
    char file[100];
    int i,j;
    double X,Y,E,tau;
    double taumax;

    snprintf(file,100,"../data/half_gaussian/tau1.txt");

    fp = fopen(file, "w");

    taumax=0.;
    for(i=0;i<L;++i){   
        X=uniform();
        Y=uniform();
        E=fabs(sqrt(-2*sigma2*log(X))*cos(2*PI*Y));
        tau=tau0*exp(E/T);
        fprintf(fp,"%lf\n",tau);
    }
  fprintf(fp,"\n sigma2=%lf",sigma2);

    fclose(fp);
}