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

int main(int argc, char* argv[]){
    std::cout<<"bloup"<<std::endl;
    KmerManipulatorACGT tmp = KmerManipulatorACGT(31);
    std::ifstream f(argv[1], std::ios::in);
    ConwayBromageBM testBM(f,&tmp);
    f.close();
    // std::ifstream f2(argv[1], std::ios::in);
    // ConwayBromageSD testSD(f2,&tmp);
    // f2.close();
    // KmerManipulatorACGT km(30);
    std::cout<<testBM.test()<<std::endl;
    // for(int i=0;i<100000;i++){
    //     if(testBM.contains(i)!=testSD.contains(i)){
    //         std::cout<<"error BM="<<testBM.contains(i)<<"and SD ="<<testSD.contains(i)<<std::endl;
    //         std::cout<<km.decode(i)<<std::endl;
    //         std::cout<<i<<std::endl;
    //     }
    // }

}
