#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <random>
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

#define L 100

int main(int argc, char *argv[]){
    FILE *fp;
    char file[100];
    int i,j,k;
    int m,m2;
    double tau[L],taumax,tauave;

    double alpha=2.5,n=4.;

    int a=atoi(argv[1]);

    tauave=pow(alpha/(alpha-1),n);//mean value of gamma distribution

    taumax=0.;
    snprintf(file,100,"../data/gamma/tau%d.txt",a);
    fp=fopen(file,"r");
    for(i=0;i<L;++i){
        fscanf(fp,"%lf\n",&tau[i]);

        if(taumax<tau[i]){
            taumax=tau[i];
            m=i;
        }
    }
    fclose(fp);

    for(i=1;i<=7;++i){
        taumax=0.;
        for(j=0;j<L;++j){
            if(j!=m){
                if(taumax<tau[j]){
                    taumax=tau[j];
                    m2=j;
                }
            }
        }

        tau[m2]=tauave;

        snprintf(file,100,"../data/gamma/tau_remove%d_%d.txt",a,i);
        fp=fopen(file,"w");
        for(j=0;j<L;++j){
            fprintf(fp,"%lf\n",tau[j]);
        }
        fclose(fp);
    }
}