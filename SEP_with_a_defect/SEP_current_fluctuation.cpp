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
#define N 10
#define M 100
#define S 1E5
#define tmax 1E5
#define taumax 10.

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

int main(){
    FILE *fp;
    char file[100];
    int i,j,k;
    int m;
    double count;
    int ds=S/100;
    double tsum,tmin,dt=tmax/M;
    double tau[L],Q[M];
    double X;

    printf("%lf %d\n",taumax,L);

    double b=0.;

    for(i=0;i<L;++i){
        tau[i]=1.;
    }

    tau[0]=taumax;

    printf("%lf\n",tau[0]);

    int x[N];
    double t[N];

    for(i=0;i<M;++i){
        Q[i]=0.;
    }

    for(i=0;i<N;++i){
        x[i]=i;
        X=uniform();
        t[i]=-tau[i]*log(X);
    }

    tsum=0.;
    while(tsum<1E3){
        tmin=t[0];
        m=0;
        for(i=1;i<N;++i){
            if(tmin>t[i]){
                tmin=t[i];
                m=i;
            }
        }

        for(i=0;i<N;++i){
            t[i]-=tmin;
        }

        tsum+=tmin;

        dynamics(x,&count,m);
        
        X=uniform();
        t[m]=-tau[x[m]]*log(X);
    }

    for(i=0;i<S;++i){
        tsum=0.;
        count=0;
        int l=0;
        while(tsum<tmax){
            tmin=t[0];
            m=0;
            for(j=1;j<N;++j){
                if(tmin>t[j]){
                    tmin=t[j];
                    m=j;
                }
            }

            for(j=0;j<N;++j){
                t[j]-=tmin;
            }

            tsum+=tmin;

            while(tsum>(l+1)*dt){
                Q[l]+=pow(count,2);
                ++l;
                if((l+1)*dt>tmax){
                    break;
                }
            }

            dynamics(x,&count,m);


            X=uniform();
            t[m]=-tau[x[m]]*log(X);
        }

        if(i%ds==0){
            snprintf(file,100,"data/current/SEP_current%d_%d.txt",(int)taumax,N);
            fp=fopen(file,"w");
            for(j=0;j<M;++j){
                fprintf(fp,"%e %lf\n",(j+1)*dt,Q[j]/(i+1));
            }
            fclose(fp);
        }
    }

    snprintf(file,100,"data/current/SEP_current%d_%d.txt",(int)taumax,N);
    fp=fopen(file,"w");
    for(i=0;i<M;++i){
        fprintf(fp,"%e %lf\n",(i+1)*dt,Q[i]/S);
    }
    fclose(fp);
}