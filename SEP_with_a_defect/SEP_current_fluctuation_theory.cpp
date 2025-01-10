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
#define taumax 100

double D(double g,double s,int N){
    double f=N*(L-N)*(-L*(s-1)*((N-1)*(s-1)+L-2)*g+s*L*((s-1)*(N-1)+L))/(L-1)/(s*N+L-N)/(-(L-2)*(s-1)*g+s*L);

    return f;
}

double g(double g1,double g2,double s,int N){
    double f=((L-4)*(L-1-N)*g1*g2-(L-3)*(L-2-N)*g1-(N-2)*g2)/((N-2)*g1+(L-3)*(L-2-N)*g2-(L-4)*(L-1-N));

    return f;
}

int main(){
    FILE *fp;
    char file[100];
    int i,j,k;
    double g1,g2;

    snprintf(file,100,"data/SEP_correlation_current_rho%d_2.txt",taumax);
    fp=fopen(file,"r");
    fscanf(fp,"%lf %lf %lf\n",&g1,&g1,&g1);
    fclose(fp);

    snprintf(file,100,"data/SEP_correlation_current_rho%d_%d.txt",taumax,L-2);
    fp=fopen(file,"r");
    printf("%s\n",file);
    fscanf(fp,"%lf %lf %lf\n",&g2,&g2,&g2);
    fclose(fp);

    g2=1.-g2;

    double a,b;

    b=taumax*(L-4)*(1.-g1)/(L-3)/(g2-g1)-taumax;
    a=g2*(b+taumax)-taumax;

    snprintf(file,100,"data/SEP_current_fluctuation_theory%d.txt",taumax);
    fp=fopen(file,"w");
    for(i=1;i<L;++i){
        fprintf(fp,"%e %e %e\n",(double)i/L,D(g(g1,g2,taumax,i),taumax,i),g(g1,g2,taumax,i));
    }
    fclose(fp);
}