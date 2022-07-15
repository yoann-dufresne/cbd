#include <ConwayBromageLib.h>
#include <sdsl/sd_vector.hpp>
#include <sdsl/vectors.hpp>
#include <bitset>
#include <queue>
#include <stack>
#include <unordered_set>
#include <list>


unordered_set<uint64_t> requestnext(int range, ConwayBromage& cb, uint64_t kmer){
    unordered_set<uint64_t> viewed;
    stack<uint64_t> kmerstack1;
    stack<uint64_t> kmerstack2;
    uint64_t mask=0;
    uint64_t kmertmp;
    //create a bit-mask to supress the first 2 bits of the kmer binary words
    for(int i=0;i<cb.getKmerManipulator()->getSize()-1;i++){
        mask=mask|(uint64_t)1<<i*2|(uint64_t)1<<((i*2)+1);
    }


    if(range){  
        uint8_t tmp=cb.successors(kmer);
        bitset<8> bloup=tmp;
        for(int i=0;i<8;i++){
            if(i<4){
                if(tmp&(1<<i)){
                    uint64_t newkmer=(kmer>>2)+((uint64_t)(3-i)<<(cb.getKmerManipulator()->getSize()-2)*2);
                    kmerstack1.push(newkmer);
                }  
            }else{
                if(tmp&(1<<i)){
                    uint64_t newkmer=((kmer<<2)&mask)+7-i;
                    kmerstack1.push(newkmer);
                }
            }
        }
    }

    for(int i=1;i<range;i++){
        std::cout<<i<<std::endl;
        if(i%2==0){
            while(!kmerstack2.empty()){
                kmertmp=kmerstack2.top();
                uint8_t succ=cb.successors(kmertmp);
                viewed.insert(kmertmp);
                for(int i=0;i<8;i++){
                    if(i<4){
                        if(succ&(1<<i)){
                            uint64_t newkmer=(kmertmp>>2)+((uint64_t)(3-i)<<(cb.getKmerManipulator()->getSize()-2)*2);
                            if (viewed.find(newkmer) == viewed.end()){//this check if the kmer has not already been seen 
                                kmerstack1.push(newkmer);

                            }
                        }  
                    }else{
                        if(succ&(1<<i)){
                            uint64_t newkmer=((kmertmp<<2)&mask)+7-i;
                            bitset<64> bloup=newkmer;
                            if(viewed.find(newkmer) == viewed.end()){//this check if the kmer has not already been seen 
                                kmerstack1.push(newkmer);
                            }                        
                        }
                    }
                }
                kmerstack2.pop();
            }
        }else{
            while(!kmerstack1.empty()){
                kmertmp=kmerstack1.top();
                uint8_t succ=cb.successors(kmertmp);
                bitset<8> bloup=succ;
                viewed.insert(kmertmp);
                for(int i=0;i<8;i++){
                    if(i<4){
                        if(succ&(1<<i)){
                            uint64_t newkmer=(kmertmp>>2)+((uint64_t)(3-i)<<(cb.getKmerManipulator()->getSize()-2)*2);
                            
                            if(viewed.find(newkmer) == viewed.end()){//this check if the kmer has not already been seen 
                                kmerstack2.push(newkmer);
                            }                        
                        }  
                    }else{
                        if(succ&(1<<i)){
                            uint64_t newkmer=((kmertmp<<2)&mask)+7-i;
                            if(viewed.find(newkmer) == viewed.end()){//this check if the kmer has not already been seen 
                                kmerstack2.push(newkmer);
                            }                        
                        }
                    }
                }
                kmerstack1.pop();

            }
 
        }
    }
    while(!kmerstack1.empty()){
        kmertmp=kmerstack1.top();
        viewed.insert(kmertmp);
        kmerstack1.pop();
    }while(!kmerstack2.empty()){
        kmertmp=kmerstack2.top();
        viewed.insert(kmertmp);
        kmerstack2.pop();
    }
    return viewed;

}
/**
 * @brief iterate over the list of kmer and call next on each of them
 */
void multnextrequest(int range,ConwayBromage& cb,list<uint64_t> kmers){
    for(uint64_t tmp : kmers){
        requestnext(range,cb,tmp);
    }
}
/*
./main file k range mult mfile/first
file the file with sorted kmer to be inserted
range the range wanted for next request
mult a boolean indicating if we want multiple request or just one, 1 mean multiple
mfile/first either the file who contain all the kmer to test next on, or a kmer in the case of only one request





*/

int main(int argc,char* argv[]){
    KmerManipulatorACGT tmp = KmerManipulatorACGT(atoi(argv[2]));
    std::ifstream f(argv[1], std::ios::in);
    ConwayBromageSD test(f,&tmp);
    KmerManipulatorACGT km(atoi(argv[2])-1);
    if(atoi(argv[4])){
        multnextrequest(atoi(argv[3]),test,)
    }else{
        requestnext(atoi(argv[3]),test,km.encode(argv[5]));
    }
    
}
