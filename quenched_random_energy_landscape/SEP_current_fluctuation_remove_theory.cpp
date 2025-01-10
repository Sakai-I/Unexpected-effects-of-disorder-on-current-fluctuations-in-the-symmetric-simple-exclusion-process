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

#define L 100 //system size

double D(double g,double s,int N){
    double f=N*(L-N)*(s*(N-1)+L+1-N-2*pow(s-1,2)*(N-1)*g/(s*L-(s-1)*(L-2)*g))/(L-1)/(s*N+L-N);

    return f;
}

double g(double g1,double g2,double s,int N){
    double f=((L-4)*(L-1-N)*g1*g2-(L-3)*(L-2-N)*g1-(N-2)*g2)/((N-2)*g1+(L-3)*(L-2-N)*g2-(L-4)*(L-1-N));

    return f;
}

//argv[1]=folder name

int main(int argc,char *argv[]){
    FILE *fp;
    char file[100];
    int i,j,k;
    int m;
    double tau[L],s,tauave,tausum;
    double g1,g2,g3;
    double epsilon;
    double A,B;

    int a=atoi(argv[2]);
    
    s=0.;
    tausum=0.;
    snprintf(file,100,"data/%s/tau%d.txt",argv[1],a);
    fp=fopen(file,"r");
    for(i=0;i<L;++i){
        fscanf(fp,"%lf\n",&tau[i]);
        tausum+=1./tau[i];
        if(s<tau[i]){
            s=tau[i];
            m=i;
        }
    }

    tauave=0.;
    for(i=0;i<L;++i){
        if(i!=m){
            tauave+=tau[i];
        }
    }

    tauave/=L-1;

    s/=tauave;

    epsilon=0.;
    for(i=0;i<L-1;++i){
        epsilon+=tau[i]*tau[i+1];
    }

    epsilon+=tau[L-1]*tau[0];

    epsilon=pow(L,2)/epsilon/tausum*tauave;

    g3=s*L*((s*(L-2)+2)/(s*(L-1)+1)-epsilon)/((s-1)*(L-2)*(s*L/(s*(L-1)+1)-epsilon));

    j=0;
    snprintf(file,100,"data/%s/SEP_with_a_defect_correlation_current_rho%d_7.txt",argv[1],a);
    fp=fopen(file,"r");
    fscanf(fp,"%lf %lf %lf\n",&A,&B,&g1);
    fscanf(fp,"%lf %lf %lf\n",&A,&B,&g2);
    fclose(fp);

    g2=1.-g2;

    tausum=0.;
    for(i=0;i<L;++i){
        tausum+=tau[i];
    }

    snprintf(file,100,"data/%s/SEP_current_fluctuation_theory_remove%d_7.txt",argv[1],a);
    fp=fopen(file,"w");
    for(i=1;i<L;++i){
        fprintf(fp,"%e %e %e\n",(double)i/L,D(g(g1,g2,s,i)*g3,s,i)/tauave,g(g1,g2,s,i)*g3);
    }
    fclose(fp);
}