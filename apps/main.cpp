#include <cstdlib>
#include <string>
#include <iostream>
#include <istream>
#include <fstream>
#include <random>
#include <chrono>
#include <bm.h>
#include <ConwayBromageLib.h>
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

int main(int argc,char* argv[]){
    KmerManipulatorACGT tmp = KmerManipulatorACGT(31);
    KmerManipulatorACGT tmp2 = KmerManipulatorACGT(30);
    std::ifstream f(argv[1], std::ios::in);
    ConwayBromageSD test(f,&tmp);
    std::cout<<test.getSequence().size()<<std::endl;
    test.serialize("./bloup");
    ConwayBromageSD test2=ConwayBromageSD::deserialize("./bloup",&tmp);
    std::cout<<test2.contains(tmp2.encode("AAAAACACTATTAGCATAAGCAGTTGTGGC"))<<std::endl;

    
}
