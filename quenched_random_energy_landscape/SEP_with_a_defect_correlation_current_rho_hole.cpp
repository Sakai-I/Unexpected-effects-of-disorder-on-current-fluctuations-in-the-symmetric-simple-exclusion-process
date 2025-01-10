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
#define tmax 1E4

void dynamics(int x[],double *Q,int N,int m,int y[]){
    int i,j,k;
    double X;

    X=uniform();
    if(y[m]==0){
        if(m==0){
            if(x[m]==0){
                if(x[N-1]!=L-1){
                    x[m]=L-1;
                    ++*Q;
                }
            }

            else if(x[m]-x[N-1]!=1){
                --x[m];
                ++*Q;
            }
        }

        else{
            if(x[m]==0){
                if(x[m-1]!=L-1){
                    x[m]=L-1;
                    ++*Q;
                }
            }

            else if(x[m]-x[m-1]!=1){
                --x[m];
                ++*Q;
            }
        }
    }

    else{
        if(m==N-1){
            if(x[m]==L-1){
                if(x[0]!=0){
                    x[m]=0;
                    --*Q;
                }
            }

            else if(x[0]-x[m]!=1){
                ++x[m];
                --*Q;
            }
        }

        else{
            if(x[m]==L-1){
                if(x[m+1]!=0){
                    x[m]=0;
                    --*Q;
                }
            }

            else if(x[m+1]-x[m]!=1){
                ++x[m];
                --*Q;
            }
        }
    }
}

int main(int argc, char *argv[]){
    FILE *fp;
    char file[100];
    int i,j,k,l,m;
    double Q;
    int ds=10000;
    double tsam,tmin;
    double taumax,tauave,s;
    double tau[L],C[L][L],tl[L],tr[L],t1,t2;
    double X;

    printf("%s %s\n",argv[1],argv[2]);

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

    tl[0]=2.*tau[L-1];
    tr[0]=2.*tau[1];
    for(i=1;i<L-1;++i){
        tl[i]=2*tau[i-1];
        tr[i]=2*tau[i+1];
    }
    tl[L-1]=2.*tau[L-2];
    tr[L-1]=2.*tau[0];

    int N=2;
    int x[N],y[N];
    double t[N];

    for(i=0;i<N;++i){
        x[i]=i;
        X=uniform();
        t1=-tl[i]*log(X);
        X=uniform();
        t2=-tr[i]*log(X);

        if(t1<t2){
            y[i]=0;//left
            t[i]=t1;
        }

        else{
            y[i]=1;
            t[i]=t2;
        }
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

        dynamics(x,&Q,N,m,y);

        X=uniform();
        t1=-tl[m]*log(X);
        X=uniform();
        t2=-tr[m]*log(X);

        if(t1<t2){
            y[m]=0;//left
            t[m]=t1;
        }

        else{
            y[m]=1;
            t[m]=t2;
        }
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

            dynamics(x,&Q,N,m,y);

            X=uniform();
            t1=-tl[x[m]]*log(X);
            X=uniform();
            t2=-tr[x[m]]*log(X);

            if(t1<t2){
                y[m]=0;//left
                t[m]=t1;
            }

            else{
                y[m]=1;
                t[m]=t2;
            }
        }

        if(i%ds==0){
            printf("%d %e %e %e\n",i/ds,(C[0][1]-C[0][L-1])/tmax/(i+1)/2,(C[1][1]-C[L-1][L-1])/tmax/(i+1)/2,(C[0][1]-C[0][L-1])/(C[1][1]-C[L-1][L-1]));
        }
    }
    printf("\n");


    printf("%s %s\n\n",argv[1],argv[2]);

    printf("%d %e %e %e\n",i/ds,(C[0][1]-C[0][L-1])/tmax/S/2,(C[1][1]-C[L-1][L-1])/tmax/S/2,(C[0][1]-C[0][L-1])/(C[1][1]-C[L-1][L-1]));
}
