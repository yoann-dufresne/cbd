#include <cstdlib>
#include <string>
#include <iostream>
#include <istream>
#include <fstream>
#include <random>
#include <chrono>
#include <ConwayBromageDatastructure.h>
#include <vector>
#include <bitset>
#include <sdsl/sd_vector.hpp>
#include <sdsl/vectors.hpp>
#include <cmath>
#include <future>


int main(int argc, char* argv[]){

    KmerManipulatorACGT tmpkm = KmerManipulatorACGT(31);
    auto a=ConwayBromageSD::deserialize(argv[1],&tmpkm);
    std::ifstream f(argv[2]);
    ConwayBromageSD b(f,&tmpkm);
    for(uint64_t i=0;i<5000000;i++){ 
        if(a.contains(i)!=b.contains(i)){
            std::cout<<"error"<<std::endl;
        }
    }
}
