#include <ConwayBromageLib.h>
#include <sdsl/sd_vector.hpp>
#include <sdsl/vectors.hpp>
#include <bitset>
#include <list>

void randomcontainsrequest(int nb,ConwayBromage& cb){
    for(int i=0;i<nb;i++){
        cb.contains(rand()%((uint64_t)(std::pow(4,(cb.getKmerManipulator()->getSize()-1) ) ) ) );
    }
}

void linearcontainsrequest(int nb,ConwayBromage& cb,int start=0,int jump=1){
    for(int i=start;i<nb+start*jump;i+=jump){
        cb.contains(i);
    }
}

void randomsuccessorsrequest(int nb,ConwayBromage& cb){
    for(int i=0;i<nb;i++){
        cb.successors(rand()%((uint64_t)(std::pow(4,(cb.getKmerManipulator()->getSize()-1) ) ) ) );
    }
}

void linearsuccessorsrequest(int nb,ConwayBromage& cb,int start=0,int jump=1){
    for(int i=start;i<nb+start*jump;i+=jump){
        cb.contains(i);
    }
}
uint64_t random64bit(){
    uint32_t km1=rand();
    uint32_t km2=rand();
    return  km1|(km2<<32);

}


void successivecontainsrequest(int nb,ConwayBromage& cb){
    uint64_t kmer=random64bit();
    uint64_t mask=0;
    for(int i=0;i<cb.getKmerManipulator()->getSize()-1;i++){
        mask=mask|(uint64_t)1<<i*2|(uint64_t)1<<i*2+1;
    }

    //choose randomly one of the next successor to 
    for(int i=0;i<nb;i++){
        cb.contains(kmer);
        kmer=(kmer<<2)&mask;
        int r=rand()%4;
        kmer+=r;
        // bitset<64> tmp(kmer);
        // std::cout<<tmp<<std::endl;
           
    }
    
}

void successivesuccessorrequest(int nb,ConwayBromage& cb){
    uint64_t kmer=random64bit();
    uint64_t mask=0;
    for(int i=0;i<cb.getKmerManipulator()->getSize()-1;i++){
        mask=mask|(uint64_t)1<<i*2|(uint64_t)1<<i*2+1;
    }

    //choose randomly one of the next successor to 
    for(int i=0;i<nb;i++){
        cb.successors(kmer);
        kmer=(kmer<<2)&mask;
        int r=rand()%4;
        kmer+=r;
    }
}

list<uint64_t> buffer(int nb,istream& f,KmerManipulator* a){
    list<uint64_t> buff;
    for(int i=0;i<nb;i++){
        std::string tmp;
        getline(f,tmp);
        buff.push_back(a->encode(tmp));
    }
    return buff;
}

void percenttest(int nb,int percent,list<uint64_t>& buffer,ConwayBromage& cb, bool contains){
    std::list<uint64_t>::iterator it =buffer.begin();
    for(int i=0;i<nb;i++){
        int p=rand()%100;
        auto a=cb.getKmerManipulator();
        std::string tmp;
        if(p>=percent){
            if(contains){
                randomcontainsrequest(1,cb);
            }else{
                randomsuccessorsrequest(1,cb);
                //std::cout<<"random"<<std::endl;

            }
        }else{
            if(contains){
                cb.contains(*it);
                it++;
            }else{
                cb.successors(*it);
                //std::cout<<a->decode(*it)<<std::endl;
                it++;
            }
        }
    }

}
void sequencetest(list<uint64_t>& buffer,ConwayBromage& cb,bool contains){
    for(uint64_t a: buffer){
        if(contains){
            cb.contains(a);
        }else{
            cb.successors(a);
        }
    }
}
void percentsequencetest(int percent,list<list<uint64_t>>& sequences,ConwayBromage& cb,bool contains){
    for(list<uint64_t> tmp:sequences){
        if(rand()%100>=percent){
            sequencetest(tmp,cb,contains);
        }else{
            if(contains){
                successivecontainsrequest(tmp.size(),cb);
            }else{
                successivesuccessorrequest(tmp.size(),cb);
            }
        }
    }

}


std::string toBinary(int n)
{
    std::string r;
    while(n!=0) {r=(n%2==0 ?"0":"1")+r; n/=2;}
    return r;
}


int main(int argc, char* argv[]){
    if((argc!=8)&&(argc!=7)&&(argc!=6)){
        std::cout<<argc<<std::endl;
        std::cout<<"./main file k nb random type percent sfile\n";
        std::cout<<"file the file wherer the kmer sorted are\n";
        std::cout<<"k the length of the kmer\n";
        std::cout<<"nb the number of the request to test\n";
        std::cout<<"random or successive kmer to test\n";
        std::cout<<"type, contain or successor\n";
        std::cout<<"percent the percentage of kmer that exist to be tested\n";
        std::cout<<"sfile the shuffled file with kmer";
        std::cout<<std::endl;
        std::ifstream f(argv[1], std::ios::in);
        KmerManipulatorACGT km(31);
        KmerManipulatorACGT km2(30);        
        ConwayBromageSD cbd(f,&km);
        std::ifstream f2(argv[2], std::ios::in);
        auto buff=buffer(200,f2,&km2);
        percenttest(200,90,buff,cbd,0);
        std::cout<<"bloup";


        exit(1);
    }
    std::ifstream f(argv[1], std::ios::in);
    KmerManipulatorACGT km(atoi(argv[2]));    
    KmerManipulatorACGT k1(atoi(argv[2])-1);
    ConwayBromageSD cbd(f,&km);
    auto start = std::chrono::steady_clock::now();
    if(argv[4]=="successive"){
        start = std::chrono::steady_clock::now();
        if(argv[5]=="contains"){
            successivecontainsrequest(atoi(argv[3]),cbd);
        }else{
            successivesuccessorrequest(atoi(argv[3]),cbd);
        }
    }else{
        if(argc!=6&&atoi(argv[6])>0){
            std::ifstream fs(argv[7]);
            auto tmp=buffer((atoi(argv[3])*2*atoi(argv[6]))/100,fs,&k1);
            start = std::chrono::steady_clock::now();
            percenttest(atoi(argv[3]),atoi(argv[6]),tmp,cbd,(argv[5]=="contains"));
        }else{
            start = std::chrono::steady_clock::now();
            if(argv[5]=="contains"){
                randomcontainsrequest(atoi(argv[3]),cbd);
            }else{
                randomsuccessorsrequest(atoi(argv[3]),cbd);
            }

        }
    }
    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;
    std::cout<<elapsed_seconds.count()<<std::endl;

    

}
