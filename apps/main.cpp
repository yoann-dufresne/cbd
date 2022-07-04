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
std::string toBinary(int n)
{
    std::string r;
    while(n!=0) {r=(n%2==0 ?"0":"1")+r; n/=2;}
    return r;
}


int main(int argc, char* argv[]){
    std::vector<bm::bvector<>> test;
    bm::bvector<> tmp;
    test.push_back(tmp);
    test[0].set(50);
    bm::bvector<> tmp2;
    test.push_back(tmp2);
    std::cout<<test[0][50]<<std::endl;
    std::cout<<test[1][50]<<std::endl;
    test.resize(50);
    Intermediate a(1);
    a.set(50);

    std::cout<<a.present(50)<<std::endl;
    KmerManipulatorACGT tmpkm = KmerManipulatorACGT(31);
    std::ifstream f(argv[1], std::ios::in);
    ConwayBromageBM test2(f,&tmpkm);
    f.close();
    KmerManipulatorACGT km(30);
    //std::cout<<toBinary(test.successors(km.encode("AAAAAGTCTGCTACTCGAAAAAAGTCTGCA")))<<std::endl;
    
}
