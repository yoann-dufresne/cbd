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


void successivecontainsrequest(int nb,ConwayBromage& cb,uint64_t kmer){
    std::cout<<"test"<<std::endl;
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
        bitset<64> tmp(kmer);
        std::cout<<tmp<<std::endl;
           
    }
    
}

void successivesuccessorrequest(int nb,ConwayBromage& cb,uint64_t kmer){
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
 /**
 * @brief 
 * 
 * @param percent the percent of random km from the file you want in your test
 * @param f a file of random kmer that we know are in cb(just shufle the original file)
 */
void percenttest(int nb,int percent,list<uint64_t> buffer,ConwayBromage& cb, bool contains){
    std::list<uint64_t>::iterator it =buffer.begin();
    for(int i=0;i<nb;i++){
        int p=rand()%100;
        auto a=cb.getKmerManipulator();
        std::string tmp;
        if(p>percent){
            if(contains){
                randomcontainsrequest(1,cb);
            }else{
                randomsuccessorsrequest(1,cb);
            }
        }else{
            if(contains){
                cb.contains(*it);
                it++;
            }else{
                cb.successors(*it);
                it++;
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
    std::string first("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");
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
        exit(1);
    }
    std::ifstream f(argv[1], std::ios::in);
    KmerManipulatorACGT km(atoi(argv[2]));    
    KmerManipulatorACGT k1(atoi(argv[2])-1);
    uint64_t first1=k1.encode(first);
    ConwayBromageSD cbd(f,&km);
    auto start = std::chrono::steady_clock::now();
    if(argv[4]=="successive"){
        auto start = std::chrono::steady_clock::now();
        if(argv[5]=="contains"){
            successivecontainsrequest(atoi(argv[3]),cbd,first1);
        }else{
            successivesuccessorrequest(atoi(argv[3]),cbd,first1);
        }
    }else{
        if(argc!=6&&atoi(argv[6])>0){
            std::ifstream fs(argv[7]);
            auto tmp=buffer((atoi(argv[3])*2*atoi(argv[6]))/100,fs,&k1);
            auto start = std::chrono::steady_clock::now();
            percenttest(atoi(argv[3]),atoi(argv[6]),tmp,cbd,(argv[5]=="contains"));
        }else{
                auto start = std::chrono::steady_clock::now();
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
