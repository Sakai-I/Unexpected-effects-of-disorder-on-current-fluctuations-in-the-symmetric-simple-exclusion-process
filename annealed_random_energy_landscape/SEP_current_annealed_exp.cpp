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
#define M 100
#define S 1E5
#define tmax 1E5

void dynamics(int x[],double t[],double T,double *tsam,double *count,int N){
    int i,j,k;
    int m;
    double X;

    double tmin=__DBL_MAX__;
    for(i=0;i<N;++i){
        if(tmin>t[i]){
            tmin=t[i];
            m=i;
        }
    }

    *tsam+=tmin;
    for(i=0;i<N;++i){
        t[i]-=tmin;
    }

    X=uniform();
    if(X<0.5){
        if(m==0){
            if(x[m]==0){
                if(x[N-1]!=L-1){
                    x[m]=L-1;
                    --*count;
                }
            }

            else if(x[m]-x[N-1]!=1){
                --x[m];
                --*count;
            }
        }

        else{
            if(x[m]==0){
                if(x[m-1]!=L-1){
                    x[m]=L-1;
                    --*count;
                }
            }

            else if(x[m]-x[m-1]!=1){
                --x[m];
                --*count;
            }
        }
    }

    else{
        if(m==N-1){
            if(x[m]==L-1){
                if(x[0]!=0){
                    x[m]=0;
                    ++*count;
                }
            }

            else if(x[0]-x[m]!=1){
                ++x[m];
                ++*count;
            }
        }

        else{
            if(x[m]==L-1){
                if(x[m+1]!=0){
                    x[m]=0;
                    ++*count;
                }
            }

            else if(x[m+1]-x[m]!=1){
                ++x[m];
                ++*count;
            }
        }
    }

    X=uniform();
    t[m]=pow(X,-1./T);
}

int main(int argc,char *argv[]){
    FILE *fp;
    char file[100];
    int i,j,k;
    int m;
    double count;
    int ds=S/100;
    double tsam,tmin,dt=tmax/M;
    double J2[M];
    double X;

    double T=atof(argv[1]); //temperature
    int N=atoi(argv[2]); //particle number
    int x[N];
    double t[N];

    snprintf(file,100,"data/current/SEP_current%.1f_%d.txt",T,N);

    for(i=0;i<M;++i){
        J2[i]=0.;
    }

    for(i=0;i<N;++i){
        x[i]=i;
        X=uniform();
        t[i]=pow(X,-1./T);
    }

    tsam=0.;
    while(tsam<1E3){
        dynamics(x,t,T,&tsam,&count,N);
    }

    for(i=0;i<S;++i){
        count=0.;
        tsam=0.;
        int l=0;
        while(tsam<tmax){
            dynamics(x,t,T,&tsam,&count,N);

            while(tsam>(l+1)*dt){
                J2[l]+=pow(count,2);
                ++l;
                if((l+1)*dt>tmax){
                    break;
                }
            }
        }

        if(i%ds==0){
            fp=fopen(file,"w");
            for(j=0;j<M;++j){
                fprintf(fp,"%e %lf\n",(j+1)*dt,J2[j]/(i+1));
            }
            fclose(fp);
        }
    }
    printf("\n");

    fp=fopen(file,"w");
    for(i=0;i<M;++i){
        fprintf(fp,"%e %lf\n",(i+1)*dt,J2[i]/S);
    }
    fclose(fp);
}