#include <ConwayBromageLib.h>
#include <sdsl/sd_vector.hpp>
#include <sdsl/vectors.hpp>
#include <bitset>
#include <queue>
#include <stack>

class kmerRange{
    public:
    uint64_t kmer;
    int currentrange;
    kmerRange(uint64_t kmer,int cr){
        kmer=kmer;
        currentrange=cr;
    }
    kmerRange(){}
};

void requestnext(int range, ConwayBromage& cb, uint64_t kmer){
    stack<kmerRange> prev;
    stack<kmerRange> next;
    uint64_t mask=0;
    //create a bit-mask to supress the first 2 bits of the kmer binary words
    for(int i=0;i<cb.getKmerManipulator()->getSize()-1;i++){
        mask=mask|(uint64_t)1<<i*2|(uint64_t)1<<i*2+1;
    }

    if(range){
        
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
    kmerRange ktmp;
    while(!prev.empty()){
        ktmp=prev.top();
        prev.pop();
        uint64_t kmer2=ktmp.kmer>>2;
        uint8_t succes=cb.successors(kmer2);
        for(int i=0;i<4;i++){
            if(succes&(1<<i)){
                if(range>ktmp.currentrange){
                    uint64_t newkmer=kmer2+(uint64_t)i<<cb.getKmerManipulator()->getSize()*2;
                    prev.push(kmerRange(newkmer,ktmp.currentrange+1));
                }   
            }

        }
    }
    while(!next.empty()){
        ktmp=next.top();
        next.pop();
        uint64_t kmer2=(ktmp.kmer)&mask;
        uint8_t succes=cb.successors(kmer2);
        for(int i=0;i<4;i++){
            if(succes&(1<<7-i)){
                if(range>ktmp.currentrange){
                    uint64_t newkmer=(kmer2<<2)+i;
                    prev.push(kmerRange(newkmer,ktmp.currentrange+1));
                }   
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
