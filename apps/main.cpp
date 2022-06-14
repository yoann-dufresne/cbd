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
    std::ofstream f2("/home/oceane/dev/test",std::ios::out|std::ios::binary);
    test.serialize(f2);
    std::ifstream f3("/home/oceane/dev/test",std::ios::in | std::ios::binary);
    ConwayBromageSD test2=ConwayBromageSD::deserialize(f3,&tmp);
    std::cout<<toBinary(test2.successors(tmp2.encode("AAAAACAAACATTGTTCGTTTAGTTGTTAA")));
    /*
    for(int i=50;i<200000;i+=9){
        std::cout<<test.successors(i)<<"=?"<<test2.successors(i)<<std::endl;
    }*/

    
}
