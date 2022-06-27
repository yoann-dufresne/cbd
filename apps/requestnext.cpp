#include <ConwayBromageLib.h>
#include <sdsl/sd_vector.hpp>
#include <sdsl/vectors.hpp>
#include <bitset>
#include <queue>

void requestnext(int range, ConwayBromageBM cb, uint64_t kmer){
    if(range){
        uint64_t mask=0;
        //create a bit-mask to supress the first 2 bits of the kmer binary words
        for(int i=0;i<cb.getKmerManipulator()->getSize()-1;i++){
            mask=mask|(uint64_t)1<<i*2|(uint64_t)1<<i*2+1;
        }
        uint8_t tmp=cb.successors(kmer);
        std::cout<<(unsigned int)tmp<<std::endl;
        uint64_t kmertmp=kmer;
        for(int i=0;i<4;i++){
            if(tmp&(1<<(7-i))){
                kmertmp=kmertmp<<2&mask;
                kmertmp+=i;
                requestnext(range-1,cb,kmertmp);
                kmertmp=kmer;
            }
        }
    }

}

int main(){
    KmerManipulatorACGT tmp = KmerManipulatorACGT(31);
    std::ifstream f("./sorted/sorted_covid", std::ios::in);
    ConwayBromageBM test(f,&tmp);
    KmerManipulatorACGT km(30);
    uint64_t tmp3=km.encode("AAAAAAAGAACAAAGACCATTGAGTACTCT");
    requestnext(50,test,tmp3);
    
}
