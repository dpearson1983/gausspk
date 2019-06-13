#include <iostream>
#include <fstream>
#include <limits>
#include <vector>
#include <string>
#include <cmath>
#include "pods.hpp"
#include "power.hpp"

power::power() {
    this->N_bins = 0;
}

power::power(int N_bins, double k_min, double k_max, double Delta_k) {
    power::init(N_bins, k_min, k_max, Delta_k);
}

void power::init(int N_bins, double k_min, double k_max, double Delta_k) {
    this->N_bins = N_bins;
    this->k_min = k_min;
    this->k_max = k_max;
    this->Delta_k = Delta_k;
    
    this->k.resize(N_bins);
    this->P.resize(N_bins);
    this->N_k.resize(N_bins);
    
    for (int i = 0; i < N_bins; ++i) {
        this->k[i] = (i + 0.5)*Delta_k + k_min;
    }
    this->first = true;
}

void power::calculate(std::vector<double2> &dk, std::vector<double> &kx, std::vector<double> &ky,
                      std::vector<double> &kz) {
    for (size_t i = 0; i < kx.size(); ++i) {
        double k_xsq = kx[i]*kx[i];
        for (size_t j = 0; j < ky.size(); ++j) {
            double k_ysq = ky[j]*ky[j];
            for (size_t l = 0; l <= kz.size()/2; ++l) {
                double k_zsq = kz[l]*kz[l];
                double k_mag = std::sqrt(k_xsq + k_ysq + k_zsq);
                
                if (k_mag >= this->k_min && k_mag < this->k_max) {
                    int index = l + (kz.size()/2 + 1)*(j + ky.size()*i);
                    int bin = (k_mag - this->k_min)/this->Delta_k;
                    this->P[bin] += dk[index].x*dk[index].x + dk[index].y*dk[index].y;
                    if (this->first) this->N_k[bin]++;
                }
            }
        }
    }
    this->first = false;
    
    for (size_t i = 0; i < this->P.size(); ++i)
        this->P[i] /= this->N_k[i];
}

void power::write(std::string file) {
    std::ofstream fout(file);
    fout.precision(std::numeric_limits<double>::digits10);
    for (size_t i = 0; i < this->P.size(); ++i) {
        fout << this->k[i] << " " << this->P[i] << "\n";
    }
    fout.close();
}
