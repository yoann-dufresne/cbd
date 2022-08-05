#include <ConwayBromageLib.h>
#include <sdsl/sd_vector.hpp>
#include <sdsl/vectors.hpp>
#include <bitset>
#include <list>
uint64_t random64bit(){
    uint32_t km1=rand();
    uint64_t km2=rand();
    return  km1|(km2<<29);

}

list<uint64_t> randomskmer(int nb){
    list<uint64_t> ret;
    for(int i=0;i<nb;i++){
        auto a=random64bit();
        ret.push_back(random64bit());
    }
    return ret;
}

void randomcontainsrequest(int nb,ConwayBromage& cb,list<uint64_t> kmer){
    auto it=kmer.begin();
    for(int i=0;i<nb;i++){
        cb.neighbours(*it);
        it++;
    }
}


void linearcontainsrequest(int nb,ConwayBromage& cb,int start=0,int jump=1){
    for(int i=start;i<nb+start*jump;i+=jump){
        cb.contains(i);
    }
}


void randomsuccessorsrequest(int nb,ConwayBromage& cb,list<uint64_t> kmer){
    auto it=kmer.begin();
    for(int i=0;i<nb;i++){
        cb.neighbours(*it);
        it++;
    }
}
void linearsuccessorsrequest(int nb,ConwayBromage& cb,int start=0,int jump=1){
    for(int i=start;i<nb+start*jump;i+=jump){
        cb.contains(i);
    }
}
list<uint64_t> randomsequenc(int size){
    list<uint64_t> ret;
    uint64_t kmer=random64bit();
    uint64_t mask=0;
    for(int i=0;i<30;i++){
        mask=mask|(uint64_t)1<<i*2|(uint64_t)1<<i*2+1;
    }
    kmer=kmer&mask;
    ret.push_back(kmer);
    //choose randomly one of the next successor to 
    for(int i=0;i<size;i++){
        kmer=(kmer<<2)&mask;
        int r=rand()%4;
        kmer+=r;
        ret.push_back(kmer);       
    }
    return ret;
}
list<list<uint64_t>> randomlistseq(int sizeseq,int nbseq){
    list<list<uint64_t>> ret;
    for(int i=0;i<nbseq;i++){
        ret.push_back(randomsequenc(sizeseq));
    }
    return ret;
}




list<uint64_t> buffer(int nb,istream& existkmer,KmerManipulator* a){
    list<uint64_t> buff;
    for(int i=0;i<nb;i++){
        std::string tmp;
        getline(existkmer,tmp);
        buff.push_back(a->encode(tmp));
    }
    return buff;
}

//do nb test with  a percentage that we know are in the structure
void percenttest(int nb,int percent,list<uint64_t>& buffer,ConwayBromage& cb, bool contains,list<uint64_t>& random){
    std::list<uint64_t>::iterator it =buffer.begin();
    auto r=random.begin();
    for(int i=0;i<nb;i++){
        int p=rand()%100;
        std::string tmp;
        if(p>=percent){
            if(contains){
                cb.contains(*r);
                r++;
            }else{
                cb.neighbours(*r);
                //std::cout<<"random"<<std::endl;
                r++;
            }
        }else{
            if(contains){
                cb.contains(*it);
                it++;
            }else{
                cb.neighbours(*it);
                //std::cout<<a->decode(*it)<<std::endl;
                it++;
            }
        }
    }

}
void sequencetest(list<uint64_t>& buffer,ConwayBromage& cb,bool contains){
    for(uint64_t a : buffer){
        if(contains){
            cb.contains(a);
        }else{
            cb.neighbours(a);
        }
    }
}
//do nb test with  a percentage that we know are in the structure
void percentsequencetest(int nb,int percent,list<list<uint64_t>>& sequences,ConwayBromage& cb,bool contains,list<list<uint64_t>>& randomseq){
    auto it=sequences.begin();
    auto it2=randomseq.begin();
    for(int i=0;i<nb;i++){
        if(rand()%100>=percent){
            sequencetest(*it2,cb,contains);
            it2++;
        }else{
            sequencetest(*it,cb,contains);
            it++;
        }
    }

}
void multiplesequencetest(list<list<uint64_t>>& buffer,ConwayBromage& cb,bool contains){
    for(auto tmp : buffer){
        sequencetest(tmp,cb,contains);
    }
}
list<list<uint64_t>> loadseqbuff(istream& f,KmerManipulator* a){
    std::string tmp;
    list<list<uint64_t>> seqbuff;
    while(getline(f,tmp)){
        seqbuff.push_back(a->encodesequence(tmp));
    }
    return seqbuff;
}

std::string toBinary(int n)
{
    std::string r;
    while(n!=0) {r=(n%2==0 ?"0":"1")+r; n/=2;}
    return r;
}


int main(int argc, char* argv[]){
    srand(time(NULL));
    if((argc!=8)&&(argc!=7)&&(argc!=6)){
        std::cout<<argc<<std::endl;
        std::cout<<"./main file k nb random type percent sfile\n";
        std::cout<<"file the file wherer the kmer sorted are\n";
        std::cout<<"k the length of the kmer\n";
        std::cout<<"nb the number of the request to test\n";
        std::cout<<"random or  kmer sequence to test\n";
        std::cout<<"type, contain or successor\n";
        std::cout<<"percent the percentage of kmer/sequence that exist to be tested\n";
        std::cout<<"sfile the shuffled file with kmer or a file containing sequence we know exist(not formally necessary but its the point of it)";
        std::cout<<std::endl;
        exit(1);
    }
    KmerManipulatorACGT km(atoi(argv[2]));    
    KmerManipulatorACGT k1(atoi(argv[2])-1);
    ConwayBromageSD cbd=ConwayBromageSD::deserialize(argv[1],&km);
    auto start = std::chrono::steady_clock::now();
    if(std::string(argv[4])=="sequence"){
        if(argc!=6){
            std::ifstream fs(argv[7]);
            auto seqbuff=loadseqbuff(fs,&k1);
            auto rs=randomlistseq((*seqbuff.begin()).size(),seqbuff.size());
            start = std::chrono::steady_clock::now();//start the chrono
            percentsequencetest(atoi(argv[3]),atoi(argv[6]),seqbuff,cbd,(std::string(argv[5])=="contains"),rs);
        }else{
            auto rs=randomlistseq(500,atoi(argv[3]));
            start = std::chrono::steady_clock::now();//start the chrono
            multiplesequencetest(rs,cbd,(std::string(argv[5])=="contains"));
        }
    }else{
        auto rk=randomskmer(atoi(argv[3]));

        if(argc!=6){
            std::ifstream fs(argv[7]);
            auto tmp=buffer(atoi(argv[3]),fs,&k1);
            start = std::chrono::steady_clock::now();
            percenttest(atoi(argv[3]),atoi(argv[6]),tmp,cbd,(std::string(argv[5])=="contains"),rk);
        }else{
            start = std::chrono::steady_clock::now();
            if(std::string(argv[5])=="contains"){
                randomcontainsrequest(atoi(argv[3]),cbd,rk);
            }else{
                randomsuccessorsrequest(atoi(argv[3]),cbd,rk);
            }

        }
    }
    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;
    std::cout<<elapsed_seconds.count()<<std::endl;

    

}
