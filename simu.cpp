// (c) Johannes Neidhart, 2020
// compile with g++ simu.cpp -lgsl -lgslcblas

#include <cmath>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_cdf.h>

#include <gsl/gsl_randist.h>
#include <unistd.h> 
#include <vector>
#include <iostream>
#include <math.h>  

// Needed for GSL RNG
long seed;
const gsl_rng_type * T;
gsl_rng * ran;


// Calculate Prbability calculates the Cumulative distribution function to have i infectuous in 
// the next time step

void calculateProbability(unsigned int S_n, double q, int I_n, std::vector < double > & P_S){
    P_S.resize(S_n);
    double Q, power;
    int coeff;
    Q = pow(q, double(I_n));
    for(unsigned int i =0; i<S_n; i++){
        P_S.at(i) =  gsl_cdf_binomial_P(double(i),(1.-Q),double(S_n));
        //std::cout<<i<<"\t"<<q<<"\t"<<1-Q<<"\t"<<S_n<<"\t"<<I_n<<"\t"<<P_S.at(i)<<std::endl;                    
    }
}

int main(){
    // setup RNG
    T = gsl_rng_default;
    ran = gsl_rng_alloc (T);
    seed = time (NULL) * getpid();
    gsl_rng_set (ran, seed);

    int S_n;            // No. of susceptible
    int startS = 150;   //No. of susceptible at start
    int I_n;            // No. of infected
    int startI = 1;     // No. of infected at start
    int n = 0;          // time
    double q = 1. - 0.05;     // Probability to not get infected

    int N = 10000;      // Statistics

    double select;      // help variable
    int h;              // help variable

    std::vector < double > P_S;     // vector storing Probability
    std::vector < double > meanRunI;// vector to store mean I
    std::vector < double > meanRunS;// vector to store mean S
    std::vector < std::vector < double > > statistics;  // array to store full runs

    statistics.resize(N);
    for(int k = 0; k<N; k++){
        statistics.at(k).push_back(startI);
        n=0;
        S_n = startS;
        I_n = startI;

        while(true){
            n++;
            calculateProbability(S_n, q, I_n, P_S);
            select = gsl_rng_uniform(ran);
            h = 0;
            for(int i = 0; i<P_S.size(); i++){       
                if(select<P_S.at(i)){
                    h = i;
                    break;
                }
            }
            I_n = h;
            S_n -= h;
            statistics.at(k).push_back(I_n);

            if(I_n == 0){
                break;
            }            
            if(S_n <= 0){
                break;
            }
        }
    }

    int nmax = 0;
    for(int i = 0; i<N; i++){
        if(statistics.at(i).size() > nmax){
            nmax = statistics.at(i).size();
        }
    }
    meanRunI.resize(nmax);
    meanRunS.resize(nmax);
    meanRunS.at(0) = startS;
    for(int j = 0; j<nmax-1; j++){
        for(int i = 0; i<N; i++){
            if((statistics.at(i).size() > j)){
                meanRunI.at(j)+=statistics.at(i).at(j);
            }else{
            }
        }
        meanRunI.at(j)/=double(N);
        if(j>0){
            meanRunS.at(j) = meanRunS.at(j-1)-meanRunI.at(j);
        }
        
        std::cout<<meanRunI.at(j)<<"\t"<<meanRunS.at(j)<<"\t"<<startS-meanRunS.at(j)<<std::endl;
    }
    return(0);
}
