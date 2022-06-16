#include <cstdlib>
#include <string>
#include <bitset>
#include <sdsl/sd_vector.hpp>
#include <sdsl/vectors.hpp>
#include <iostream>
#include <istream>
#include <ostream>

#include "ConwayBromageLib.h"
#include "lest.hpp"

const lest::test module[] =
{
    CASE("serialize and deserialize from a ConwayBromage CBD" "[serialize]"){
        std::ifstream f("/home/oceane/dev/sortedcovid31", ios::in);
        KmerManipulatorACGT k(31);
        KmerManipulatorACGT k3(30);
        ConwayBromageBM cb(f, &k);

        std::ofstream o("/home/oceane/dev/test3",std::ios::out);
        cb.serialize(o);
        std::ifstream f2("/home/oceane/dev/test3", ios::in);
        ConwayBromageBM cb2=ConwayBromageBM::deserialize(f2, &k);
        for(int i =1  ; i < 50000000 ; i++){
                bitset<8> bitForm((unsigned)cb.successors(i));  //uint8_t of successors, bit version
                bitset<8> bitForm2((unsigned)cb2.successors(i));
                if(bitForm!=0){
                    std::cout<<bitForm<<"=?"<<bitForm2<<std::endl;
                }
                EXPECT(bitForm==bitForm2);
            }
        
    }
};

extern lest::tests & specification();

MODULE( specification(), module )
