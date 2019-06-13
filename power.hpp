#ifndef _POWER_HPP_
#define _POWER_HPP_

#include <vector>
#include "pods.hpp"

class power{
    std::vector<double> k, P;
    std::vector<long> N_k;
    int N_bins;
    double k_min, k_max, Delta_k;
    bool first;
    
    public:
        power();
        
        power(int N_bins, double k_min, double k_max, double Delta_k);
        
        void init(int N_bins, double k_min, double k_max, double Delta_k);
        
        void calculate(std::vector<double2> &dk, std::vector<double> &kx, std::vector<double> &ky,
                  std::vector<double> &kz);
        
        void write(std::string file);
        
};

#endif
