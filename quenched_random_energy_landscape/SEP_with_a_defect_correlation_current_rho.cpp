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

uint64_t get_rand(){
    return mt64();
}

#define L 100
#define S 1E7
#define N 2
#define tmax 1E4

void dynamics(int x[],double *Q,int m){
    int i,j,k;
    double X;

    X=uniform();
    if(X<0.5){
        if(m==0){
            if(x[m]==0){
                if(x[N-1]!=L-1){
                    x[m]=L-1;
                    --*Q;
                }
            }

            else if(x[m]-x[N-1]!=1){
                --x[m];
                --*Q;
            }
        }

        else{
            if(x[m]==0){
                if(x[m-1]!=L-1){
                    x[m]=L-1;
                    --*Q;
                }
            }

            else if(x[m]-x[m-1]!=1){
                --x[m];
                --*Q;
            }
        }
    }

    else{
        if(m==N-1){
            if(x[m]==L-1){
                if(x[0]!=0){
                    x[m]=0;
                    ++*Q;
                }
            }

            else if(x[0]-x[m]!=1){
                ++x[m];
                ++*Q;
            }
        }

        else{
            if(x[m]==L-1){
                if(x[m+1]!=0){
                    x[m]=0;
                    ++*Q;
                }
            }

            else if(x[m+1]-x[m]!=1){
                ++x[m];
                ++*Q;
            }
        }
    }
}

int main(int argc,char *argv[]){
    FILE *fp;
    char file[100];
    int i,j,k,l,m;
    int a;
    double Q;
    int ds=10000;
    double tsam,tmin;
    double taumax,tauave,s;
    double tau[L],C[L][L];
    double X;

    snprintf(file,100,"data/%s/tau%s.txt",argv[1],argv[2]);
    fp=fopen(file,"r");
    for(i=0;i<L;++i){
        fscanf(fp,"%lf\n",&tau[i]);
    }
    fclose(fp);

    taumax=tau[0];
    j=0;
    for(i=1;i<L;++i){
        if(taumax<tau[i]){
            taumax=tau[i];
            j=i;
        }
    }

    tauave=0.;
    for(i=0;i<L;++i){
        if(i!=j){
            tauave+=tau[i];
        }
    }

    tauave/=L-1;

    s=taumax/tauave;

    printf("%e\n",s);

    for(i=0;i<L;++i){
        tau[i]=1.;
    }
    tau[0]=s;

    int x[N];
    double t[N];

    for(i=0;i<N;++i){
        x[i]=i;
        X=uniform();
        t[i]=-tau[i]*log(X);
    }

    tsam=0.;
    while(tsam<1E4){
        double tmin=__DBL_MAX__;
        for(i=0;i<N;++i){
            if(tmin>t[i]){
                tmin=t[i];
                m=i;
            }
        }

        tsam+=tmin;
        for(i=0;i<N;++i){
            t[i]-=tmin;
        }

        dynamics(x,&Q,m);

        X=uniform();
        t[m]=-tau[x[m]]*log(X);
    }

    for(i=0;i<L;++i){
        for(j=0;j<L;++j){
            C[j][i]=0.;
        }
    }

    for(i=0;i<S;++i){
        tsam=0.;
        Q=0.;
        while(tsam<tmax){
            double tmin=__DBL_MAX__;
            for(j=0;j<N;++j){
                if(tmin>t[j]){
                    tmin=t[j];
                    m=j;
                }
            }

            tsam+=tmin;
            for(j=0;j<N;++j){
                t[j]-=tmin;
            }

            for(j=0;j<N;++j){
                for(k=0;k<N;++k){
                    C[x[j]][x[k]]+=tmin*Q;
                }
            }

            if(tsam>tmax){
                break;
            }

            dynamics(x,&Q,m);

            X=uniform();
            t[m]=-tau[x[m]]*log(X);
        }

        if(i%ds==0){
            printf("%d %e %e %e\n",i/ds,(C[0][1]-C[0][L-1])/tmax/(i+1)/2,(C[1][1]-C[L-1][L-1])/tmax/(i+1)/2,(C[0][1]-C[0][L-1])/(C[1][1]-C[L-1][L-1]));
        }
    }
    printf("\n");

    printf("%s %s\n\n",argv[1],argv[2]);

    printf("%e %e %e %e\n",tmax,(C[0][1]-C[0][L-1])/tmax/S/2,(C[1][1]-C[L-1][L-1])/tmax/S/2,(C[0][1]-C[0][L-1])/(C[1][1]-C[L-1][L-1]));
}
