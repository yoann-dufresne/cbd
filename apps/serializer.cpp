#include <cstdlib>
#include <string>
#include <iostream>
#include <istream>
#include <fstream>
#include <random>
#include <chrono>
#include <bm64.h>
#include <ConwayBromageLib.h>
#include <vector>
#include <bitset>
#include <sdsl/sd_vector.hpp>
#include <sdsl/vectors.hpp>
#include <cmath>
#include <future>


int main(int argc, char* argv[]){
    KmerManipulatorACGT tmpkm = KmerManipulatorACGT(31);
    std::ifstream f(argv[1], std::ios::in);
    ConwayBromageSD test2(f,&tmpkm);
    f.close();
    test2.serialize(argv[2]);
}
