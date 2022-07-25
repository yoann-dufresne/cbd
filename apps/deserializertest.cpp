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
    auto a=ConwayBromageSD::deserialize(argv[1],&tmpkm);
}
