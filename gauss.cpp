#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>
#include <limits>
#include <gsl/gsl_spline.h>
#include "pods.hpp"
#include "gauss.hpp"

std::vector<double> gausspk::fftFrequencies(int n, double l) {
    std::vector<double> k(n);
    double deltak = 2.0*M_PI/l;
    for (int i = 0; i <= n/2; ++i)
        k[i] = i*deltak;
    for (int i = n/2 + 1; i < n; ++i)
        k[i] = (i - n)*deltak;
    return k;
}

gausspk::gausspk() : gen((std::random_device())()), dist(0.0, 1.0) {
    this->N = {0, 0, 0, 0};
}

gausspk::~gausspk() {
    gsl_spline_free(this->Pk);
    gsl_interp_accel_free(this->acc);
}

gausspk::gausspk(int N, double L) : gen((std::random_device())()), dist(0.0, 1.0) {
    gausspk::init(N, L);
}

gausspk::gausspk(int3 N, double3 L) : gen((std::random_device())()), dist(0.0, 1.0) {
    gausspk::init(N, L);
}

gausspk::gausspk(int N, double L, std::string pkFile) : gen((std::random_device())()), dist(0.0, 1.0) {
    gausspk::init(N, L, pkFile);
}

gausspk::gausspk(int3 N, double3 L, std::string pkFile) : gen((std::random_device())()), dist(0.0, 1.0) {
    gausspk::init(N, L, pkFile);
}

void gausspk::init(int N, double L) {
    this->N = {N, N, N, N*N*N};
    this->L = {L, L, L, L*L*L};
    
    this->dk.resize(N*N*(N/2 + 1));
    
    this->kx = gausspk::fftFrequencies(N, L);
    this->ky = gausspk::fftFrequencies(N, L);
    this->kz = gausspk::fftFrequencies(N, L);
}

void gausspk::init(int3 N, double3 L) {
    this->N = {N.x, N.y, N.z, N.x*N.y*N.z};
    this->L = {L.x, L.y, L.z, L.x*L.y*L.z};
    
    this->dk.resize(N.x*N.y*(N.z/2 + 1));
    
    this->kx = gausspk::fftFrequencies(N.x, L.x);
    this->ky = gausspk::fftFrequencies(N.y, L.y);
    this->kz = gausspk::fftFrequencies(N.z, L.z);
}

void gausspk::init(int N, double L, std::string pkFile) {
    this->N = {N, N, N, N*N*N};
    this->L = {L, L, L, L*L*L};
    
    this->dk.resize(N*N*(N/2 + 1));
    
    this->kx = gausspk::fftFrequencies(N, L);
    this->ky = gausspk::fftFrequencies(N, L);
    this->kz = gausspk::fftFrequencies(N, L);
    
    if (std::ifstream(pkFile)) {
        std::vector<double> k, P;
        std::ifstream fin(pkFile);
        while (!fin.eof()) {
            double kt, Pt;
            fin >> kt >> Pt;
            if (!fin.eof()) {
                k.push_back(kt);
                P.push_back(Pt);
            }
        }
        fin.close();
        
        this->Pk = gsl_spline_alloc(gsl_interp_cspline, P.size());
        this->acc = gsl_interp_accel_alloc();
        
        gsl_spline_init(this->Pk, k.data(), P.data(), P.size());
    } else {
        std::stringstream errMsg;
        errMsg << "Cannot open input power spectrum file: " << pkFile << "\n";
        throw std::runtime_error(errMsg.str());
    }
}

void gausspk::init(int3 N, double3 L, std::string pkFile) {
    this->N = {N.x, N.y, N.z, N.x*N.y*N.z};
    this->L = {L.x, L.y, L.z, L.x*L.y*L.z};
    
    this->dk.resize(N.x*N.y*(N.z/2 + 1));
    
    this->kx = gausspk::fftFrequencies(N.x, L.x);
    this->ky = gausspk::fftFrequencies(N.y, L.y);
    this->kz = gausspk::fftFrequencies(N.z, L.z);
    
    if (std::ifstream(pkFile)) {
        std::vector<double> k, P;
        std::ifstream fin(pkFile);
        while (!fin.eof()) {
            double kt, Pt;
            fin >> kt >> Pt;
            if (!fin.eof()) {
                k.push_back(kt);
                P.push_back(Pt);
            }
        }
        fin.close();
        
        this->Pk = gsl_spline_alloc(gsl_interp_cspline, P.size());
        this->acc = gsl_interp_accel_alloc();
        
        gsl_spline_init(this->Pk, k.data(), P.data(), P.size());
    } else {
        std::stringstream errMsg;
        errMsg << "Cannot open input power spectrum file: " << pkFile << "\n";
        throw std::runtime_error(errMsg.str());
    }
}

void gausspk::setPk(std::string pkFile) {
    if (std::ifstream(pkFile)) {
        std::vector<double> k, P;
        std::ifstream fin(pkFile);
        while (!fin.eof()) {
            double kt, Pt;
            fin >> kt >> Pt;
            if (!fin.eof()) {
                k.push_back(kt);
                P.push_back(Pt);
            }
        }
        fin.close();
        
        this->Pk = gsl_spline_alloc(gsl_interp_cspline, P.size());
        this->acc = gsl_interp_accel_alloc();
        
        gsl_spline_init(this->Pk, k.data(), P.data(), P.size());
    } else {
        std::stringstream errMsg;
        errMsg << "Cannot open input power spectrum file: " << pkFile << "\n";
        throw std::runtime_error(errMsg.str());
    }
}

void gausspk::setPk(std::vector<double> &k, std::vector<double> &P) {
    this->Pk = gsl_spline_alloc(gsl_interp_cspline, P.size());
    this->acc = gsl_interp_accel_alloc();
    
    gsl_spline_init(this->Pk, k.data(), P.data(), P.size());
}

void gausspk::sample() {
    for (int i = 0; i < this->N.x; ++i) {
        int i2 = (2*this->N.x - i) % this->N.x;
        double k_xsq = this->kx[i]*this->kx[i];
        for (int j = 0; j < this->N.x; ++j) {
            int j2 = (2*this->N.y - j) % this->N.y;
            double k_ysq = this->ky[j]*this->ky[j];
            for (int k = 0; k < this->N.z/2; ++k) {
                double k_mag = std::sqrt(k_xsq + k_ysq + this->kz[k]*this->kz[k]);
                int index = k + (this->N.z/2 + 1)*(j + this->N.y*i);
                
                if (k_mag > 0) {
                    double P = gsl_spline_eval(this->Pk, k_mag, this->acc);
                    if ((i == 0 || i == N.x/2) && (j == 0 || j == N.y/2) && (k == 0 && k == N.z/2)) {
                        dk[index].x = sqrt(P)*this->dist(this->gen);
                    } else if (k == 0 || k == N.z/2) {
                        int index2 = k + (this->N.z/2 + 1)*(j2 + this->N.y*i2);
                        dk[index].x = sqrt(P/2.0)*this->dist(this->gen);
                        dk[index].y = sqrt(P/2.0)*this->dist(this->gen);
                        
                        dk[index2].x = dk[index].x;
                        dk[index2].y = dk[index].y;
                    } else {
                        dk[index].x = sqrt(P/2.0)*this->dist(this->gen);
                        dk[index].y = sqrt(P/2.0)*this->dist(this->gen);
                    }
                } else {
                    dk[index].x = 0.0;
                    dk[index].y = 0.0;
                }
            }
        }
    }
}

void gausspk::writeTXT(std::string file) {
    std::ofstream fout(file);
    fout.precision(std::numeric_limits<double>::digits10);
    for (int i = 0; i < this->N.x; ++i) {
        for (int j = 0; j < this->N.y; ++j) {
            for (int k = 0; k < this->N.z; ++k) {
                int index = k + (this->N.z/2 + 1)*(j + this->N.y*i);
                fout << i << " " << j << " " << k << " " << this->dk[index].x << " " << this->dk[index].y << "\n";
            }
        }
    }
    fout.close();
}

void gausspk::writeBIN(std::string file) {
    std::ofstream fout(file, std::ios::out|std::ios::binary);
    fout.write((char *)this->dk.data(), this->dk.size()*sizeof(double2));
    fout.close();
}

std::vector<double2> gausspk::getdk() {
    return this->dk;
}

void gausspk::getFrequencies(std::vector<double> &kx, std::vector<double> &ky, std::vector<double> &kz) {
    kx = this->kx;
    ky = this->ky;
    kz = this->kz;
}
