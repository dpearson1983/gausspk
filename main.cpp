#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include "pods.hpp"
#include "gauss.hpp"
#include "power.hpp"
#include "harppi.hpp"

std::string filename(std::string base, int digits, int num, std::string ext) {
    std::stringstream file;
    file << base << std::setw(digits) << std::setfill('0') << num << "." << ext;
    return file.str();
}

int main(int argc, char *argv[]) {
    parameters p(argv[1]);
    p.print();
    
    int3 N = {p.geti("Nx"), p.geti("Ny"), p.geti("Nz")};
    double3 L = {p.getd("Lx"), p.getd("Ly"), p.getd("Lz")};
    
    gausspk delta(N, L, p.gets("pkFile"));
    power Pk(p.geti("N_bins"), p.getd("k_min"), p.getd("k_max"), p.getd("Delta_k"));
    
    for (int real = p.geti("startNum"); real < p.geti("startNum") + p.geti("numReals"); ++real) {
        std::cout << "Realization #" << real << std::endl;
        std::string outFile = filename(p.gets("outbase"), p.geti("digits"), real, p.gets("outext"));
        
        delta.sample();
        Pk.calculate(delta.dk, delta.kx, delta.ky, delta.kz);
        Pk.write(outFile);
    }
    
    return 0;
}
