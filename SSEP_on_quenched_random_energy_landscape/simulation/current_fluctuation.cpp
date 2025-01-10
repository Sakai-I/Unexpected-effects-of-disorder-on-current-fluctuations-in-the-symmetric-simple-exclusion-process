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

#define L 100 //system size
#define N 10 //particle number
#define M 100
#define S 1E5
#define tmax 1E5 //observation time

void dynamics(int x[],double t[],double tau[],double *tsam,double *count){//dynamics
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
    t[m]=-tau[x[m]]*log(X);
}

//argv[1]=file name

int main(int argc,char *argv[]){
    FILE *fp;
    char file[100];
    int i,j,k;
    int m;
    double count;
    int ds=S/100;
    int x[N];
    double t[N];
    double tsum,tmin,dt=tmax/M;
    double tau[L],Q2[M];
    double X;

    int a=atoi(argv[2]);

    snprintf(file,100,"%s/tau%d.txt",argv[1],a);
    fp=fopen(file,"r");
    for(i=0;i<L;++i){
        fscanf(fp,"%lf\n",&tau[i]);
    }
    fclose(fp);

    for(i=0;i<M;++i){
        Q2[i]=0.;
    }

    for(i=0;i<N;++i){
        x[i]=i;
        X=uniform();
        t[i]=-tau[i]*log(X);
    }

    for(i=0;i<S;++i){
        for(j=0;j<N;++j){
            X=uniform();
            t[j]=-tau[x[j]]*log(X);
        }

        tsum=0.;
        while(tsum<1E3){//initial condition
            dynamics(x,t,tau,&tsum,&count);
        }

        count=0.;

        tsum=0.;
        int l=0;
        while(tsum<tmax){
            dynamics(x,t,tau,&tsum,&count);

            while(tsum>(l+1)*dt){
                Q2[l]+=pow(count,2);
                ++l;
                if((l+1)*dt>tmax) break;
            }
        }

        if(i%ds==0){
            printf("\b");
            printf("=>");
            fflush(stdout);

            snprintf(file,100,"%s/current/SEP_current%d_%d.txt",argv[1],a,N);
            fp=fopen(file,"w");
            for(j=0;j<M;++j){
                fprintf(fp,"%e %lf\n",(j+1)*dt,Q2[j]/(i+1));
            }
            fclose(fp);
        }
    }
    printf("\n");

    snprintf(file,100,"%s/current/SEP_current%d_%d.txt",argv[1],a,N);
    fp=fopen(file,"w");
    for(i=0;i<M;++i){
        fprintf(fp,"%e %lf\n",(i+1)*dt,Q2[i]/S);
    }
    fclose(fp);
}
