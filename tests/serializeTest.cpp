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
        std::ifstream f("./sortACGT.txt", ios::in);
        KmerManipulatorACGT k(4);
        KmerManipulatorACGT k3(3);
        ConwayBromageBM cb(f, &k);

        std::ofstream o("/home/oceane/dev/test3",std::ios::out|std::ios::binary);
        cb.serialize(o);
        std::ifstream f2("/home/oceane/dev/test3", ios::in);
        ConwayBromageBM cb2=ConwayBromageBM::deserialize(f2, &k);
        string call[4]{"A", "C", "G", "T"}; //build comrades manually
        for(int i = 1 ; i < 64 ; i++){
                bitset<8> bitForm((unsigned)cb.successors(i));  //uint8_t of successors, bit version
                bitset<8> bitForm2((unsigned)cb2.successors(i));
                std::cout<<(bitForm==bitForm2)<<std::endl;
            }
        bitset<8> bitForm((unsigned)cb.successors(1));  //uint8_t of successors, bit version
        bitset<8> bitForm2((unsigned)cb2.successors(1));
        EXPECT(bitForm==bitForm2);
    }
};

extern lest::tests & specification();

MODULE( specification(), module )
