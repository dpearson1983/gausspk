#ifndef _GAUSS_HPP_
#define _GAUSS_HPP_

#include <vector>
#include <string>
#include <random>
#include <gsl/gsl_spline.h>
#include "pods.hpp"

class gausspk{
    int4 N;
    double4 L;
    gsl_spline *Pk;
    gsl_interp_accel *acc;
    std::mt19937_64 gen;
    std::normal_distribution<double> dist;
    
    std::vector<double> fftFrequencies(int n, double l);
    
    public:
        // Some public data members for ease of use
        std::vector<double2> dk;
        std::vector<double> kx, ky, kz;
        
        gausspk();
        
        ~gausspk();
        
        gausspk(int N, double L);
        
        gausspk(int3 N, double3 L);
        
        gausspk(int N, double L, std::string pkFile);
        
        gausspk(int3 N, double3 L, std::string pkFile);
        
        void init(int N, double L);
        
        void init(int3 N, double3 L);
        
        void init(int N, double L, std::string pkFile);
        
        void init(int3 N, double3 L, std::string pkFile);
        
        void setPk(std::string pkFile);
        
        void setPk(std::vector<double> &k, std::vector<double> &P);
        
        void sample();
        
        void writeTXT(std::string file);
        
        void writeBIN(std::string file);
        
        std::vector<double2> getdk();
        
        void getFrequencies(std::vector<double> &kx, std::vector<double> &ky, std::vector<double> &kz);
        
};

#endif
