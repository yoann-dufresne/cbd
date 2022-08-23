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
std::string toBinary(int n)
{
    std::string r;
    while(n!=0) {r=(n%2==0 ?"0":"1")+r; n/=2;}
    return r;
}


int main(int argc, char* argv[]){
    KmerManipulatorACGT tmpkm = KmerManipulatorACGT(31);
    std::ifstream f(argv[1], std::ios::in);
    ConwayBromageSD test2(f,&tmpkm);
    f.close();
    KmerManipulatorACGT km(30);
    //std::cout<<toBinary(test.successors(km.encode("AAAAAGTCTGCTACTCGAAAAAAGTCTGCA")))<<std::endl;
    
}
