#include <cstdlib>
#include <string>
#include <iostream>
#include <istream>
#include <fstream>
#include <random>
#include <chrono>
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

int main(){
    KmerManipulatorACTG tmp = KmerManipulatorACTG(31);
    std::ifstream f("./test2", std::ios::in);
    ConwayBromageSD test(f,&tmp);
    KmerManipulatorACTG km(30);
    uint64_t tmp3=km.encode("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");
    uint8_t test2=test.successors(tmp3);
    printf("%x\n",test2);
    std::string suc=toBinary(test2);
    std::cout<< suc<<std::endl;
    std::cout<<std::hex<<(unsigned int)test2<<std::endl;
    std::cout<<test.contains(km.getCanonical(km.encode("CACCTACCTGGCGATTATGCGCGGTTACGC")))<<std::endl;

    //std::cout<<tmp.decode(tmp.getCanonical(tmp.encode("GAAA")))<<std::endl;

}
