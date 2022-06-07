//
// Created by Alexandra on 18/06/2020.
//

/* Transform sequences which are contain in a file in a sd_vector
 * Merge the reverse complement, fastest ASCII version for reverse complement and encode/decode
 * @param path - a string which is the path to the file which contains the generated sequence
 * @return a sd_vector which contains an encoding merging version of the sequence of the file and of the reverse complement of the sequence
 */
sd_vector<>fromFileToSdVectorWithReverseEcoli(string path){
    ifstream file(path, ios::in);  // Reading of the file which contains k-mers sequences
    if(file){   // File is now open
        string word;
        string line("");
        file >> word;   //Take the first word to analyze size of one k-mer
        file.seekg(0, ios::beg);    //Return to the beginning of the file
        uint64_t myWordLen(word.size()); //Size of k_mer, it is the 'k'
        uint64_t myOneLen(0);
        while(getline(file, line)){ //Counts the number of ones in the file
            myOneLen++;
        }
        file.clear();
        file.seekg(0, ios::beg);    //Return to the beginning of the file
        cout << "length of ones : " << myOneLen << endl;
        cout << "length of a seq : " << myWordLen << endl;
        uint64_t myTotalLen(pow(ALPHABET,myWordLen));   //Creation of the total length to create sd_vector_builders
        int_vector<> reverser(myOneLen, 0); //to simplify the reverse complement build                                                              //SLOW DOWN
        cout << "Total length : " << myTotalLen << endl;
        int i = 0;
        sd_vector_builder constructSparse(myTotalLen, myOneLen);    //A size of myTotalLen, contains myOneLen ones, right version
        sd_vector_builder constructReverse(myTotalLen, myOneLen);   //same total size and ones size, reverse version
        while(file >> word){
            constructSparse.set(encodeEcoli(word, myWordLen)); //filled to one each element which is represent by the encoding version of the sequence
            reverser[i] = reverseComplementGATBLibEcoli(encodeEcoli(word, myWordLen), myWordLen);  //each reverses of each words, compressed version                                           //SLOW DOWN
            i++;
            file >> word;
        }
        sort(reverser.begin(), reverser.end()); //sorting of the reverse otherwise it will be impossible to use sd_vector_builder                               //SLOW DOWN
        for(int i = 0 ; i < myOneLen ; i++){
            constructReverse.set(reverser[i]);  //filled for the reverse version                                                                                //SLOW DOWN
        }
        vector<sd_vector<>> seq;
        sd_vector<>finalSparseRight(constructSparse);    //Construction of the final sd_vector for the right sequence
        sd_vector<>finalSparseReverse(constructReverse);   //Construction of the final sd_vector for the reverse sequence
        sd_vector<>finalAll = merge(finalSparseRight, finalSparseReverse, finalSparseRight.size(), finalSparseReverse.size());
        file.close();
        return finalAll;
    }else{
        cout << "Error while opening" << endl;
    }
    return bit_vector{0};
}

/* Transform sequences which are contain in a file in a sd_vector
 * Merge with the reverse complement
 * @param path - a string which is the path to the file which contains the generated sequence
 * @return a sd_vector which contains an encoding merging version of the sequence of the file and of the reverse complement of the sequence
 */
sd_vector<>fromFileToSdVectorWithReverse(string path){
    ifstream file(path, ios::in);  // Reading of the file which contains k-mers sequences
    if(file){   // File is now open
        string word;
        string line("");
        file >> word;   //Take the first word to analyze size of one k-mer
        file.seekg(0, ios::beg);    //Return to the beginning of the file
        uint64_t myWordLen(word.size()); //Size of k_mer, it is the 'k'
        uint64_t myOneLen(0);
        while(getline(file, line)){ //Counts the number of ones in the file
            myOneLen++;
        }
        file.clear();
        file.seekg(0, ios::beg);    //Return to the beginning of the file
        cout << "length of ones : " << myOneLen << endl;
        cout << "length of a seq : " << myWordLen << endl;
        uint64_t myTotalLen(pow(ALPHABET,myWordLen));   //Creation of the total length to create sd_vector_builders
        int_vector<> reverser(myOneLen, 0); //to simplify the reverse complement build
        cout << "Total length : " << myTotalLen << endl;
        int i = 0;
        sd_vector_builder constructSparse(myTotalLen, myOneLen);    //A size of myTotalLen, contains myOneLen ones, right version
        sd_vector_builder constructReverse(myTotalLen, myOneLen);   //same total size and ones size, reverse version
        while(file >> word){
            if(word != "1"){
                constructSparse.set(encode(word, myWordLen)); //filled to one each element which is represent by the encoding version of the sequence
                reverser[i] = reverseComplementLexico(encode(word, myWordLen), myWordLen);  //each reverses of each words, compressed version
                i++;
            }
        }
        sort(reverser.begin(), reverser.end()); //sorting of the reverse otherwise it will be impossible to use sd_vector_builder
        for(int i = 0 ; i < myOneLen ; i++){
            constructReverse.set(reverser[i]);  //filled for the reverse version
        }
        vector<sd_vector<>> seq;
        sd_vector<>finalSparseRight(constructSparse);    //Construction of the final sd_vector for the right sequence
        sd_vector<>finalSparseReverse(constructReverse);   //Construction of the final sd_vector for the reverse sequence
        sd_vector<>finalAll = merge(finalSparseRight, finalSparseReverse, finalSparseRight.size(), finalSparseReverse.size());  //merging
        file.close();
        return finalAll;
    }else{
        cout << "Error while opening" << endl;
    }
    return bit_vector{0};
}

/* Transform sequences which are contain in a file in a sd_vector
 * Need a string which is the path to the file
 * Return a sd_vector which contains encoding version of sequences of the file
 */
sd_vector<>fromFileToSdVectorEcoli(string path){
    ifstream file(path, ios::in);  // Reading of the file which contains k-mers sequences
    if(file){   // File is now open
        string word;
        string line("");
        file >> word;   //Take the first word to analyze size of one k-mer
        file.seekg(0, ios::beg);    //Return to the beginning of the file
        int myWordLen(word.size()); //Size of k_mer, it is the 'k'
        int myOneLen(0);
        while(getline(file, line)){ //Counts the number of ones in the file
            myOneLen++;
        }
        file.clear();
        file.seekg(0, ios::beg);    //Return to he beginning of the file
        cout << "length of ones : " << myOneLen << endl;
        cout << "length of a seq : " << myWordLen << endl;
        uint64_t myTotalLen(pow(ALPHABET,myWordLen));   //Creation of the total length to create the sd_vector_builder
        cout << "Total length : " << myTotalLen << endl;
        sd_vector_builder constructSparse(myTotalLen, myOneLen);    //A size of myTotalLen, contains myOneLen ones
        while(file >> word){
            //cout << "The seq is : " << word << endl;
            constructSparse.set(encodeEcoli(word, myWordLen)); //filled to one each element which is represent by the encoding version of the sequence
            file >> word;
        }
        sd_vector<>finalSparse(constructSparse);    //Construction of the final sd_vector
        file.close();
        return finalSparse;
    }else{
        cout << "Error while opening" << endl;
    }
    return bit_vector{0};
}

/**
 * Transform a file which contains k-mers in an sd_vector.
 * This version, generates a txt file which contains the values of reverse complement.
 * @param path of the file
 * @return an sd_vector of size 4^(P-1)
 */
sd_vector<> fromFileToSdVector_TXTversion(string path){
    ifstream file(path, ios::in);  // Reading of the file which contains k-mers sequences
    if(!file) {
        cout << "Error while opening" << endl;
        return bit_vector{0};
    }

    //obtention of the k-mers' size
    string kmer("");
    file >> kmer;
    file.seekg(0, ios::beg);
    int K = kmer.size();

    //counting the number of k-mer in the file
    uint64_t nb_of_kmer = 0;
    string line("");
    while(getline(file, line)) nb_of_kmer++;
    cout << "nb_of_kmer : " << nb_of_kmer << endl;
    file.clear();
    file.seekg(0, ios::beg);
    sd_vector_builder sdv_builder_classic(pow(ALPHABET,K), nb_of_kmer);

    //first parsing to get the classic read
    while(file >> kmer){
        sdv_builder_classic.set(encode(kmer, K));//reverseComplementLexico(i_classic, K);
        file >> kmer;
    }
    file.clear();
    file.seekg(0, ios::beg);
    sd_vector<> classic(sdv_builder_classic);
    cout << "Length of classic : " << classic.size() << endl;
    cout << "size of classic (MB) : " << size_in_mega_bytes(classic) << endl;

    cout << "Completion of distinct_rev_comp.txt." << endl;
    auto b1 = std::chrono::high_resolution_clock::now();
    //second parsing to get the reverse complement
    ofstream output;
    output.open("distinct_rev_comp.txt");
    uint64_t nb_of_distinct_rev_comp = 0;
    while(file >> kmer){
        uint64_t i_classic = encode(kmer, K);
        uint64_t i_rc = reverseComplementLexico(i_classic, K);
        if(classic[i_rc]==0){ //if the reverse complement is not present in the classical reads then add it in the output
            output << i_rc;
            output << '\n';
            nb_of_distinct_rev_comp++;
        }
        file >> kmer;
    }
    output.close();
    auto e1 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed1 = e1 - b1;
    cout << "--> Time (s): " << elapsed1.count() << endl;

    //sort the output
    cout << "distinct_rev_comp.txt completed. Now we sort." << endl;
    auto b2 = std::chrono::high_resolution_clock::now();
    string command = "sort -n distinct_rev_comp.txt -o distinct_rev_comp.txt";
    system(command.c_str());
    cout << "Sort finished. " << endl;
    auto e2 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed2 = e2 - b2;
    cout << "--> Time (s): " << elapsed2.count() << endl;
    cout << "Now we build the final sd_vector." << endl;

    //now we build the sd_vector
    auto b3 = std::chrono::high_resolution_clock::now();
    ifstream input;
    input.open("distinct_rev_comp.txt");  // Reading of the file which contains k-mers sequences
    if(!input) {
        cout << "Error while opening distinct_rev_comp.txt" << endl;
        return bit_vector{0};
    };
    sd_vector_builder sdvB(pow(ALPHABET,K), nb_of_kmer + nb_of_distinct_rev_comp);
    uint64_t i = 0, i_classic = 1, i_rc = 0, last_rc = 0;

    sd_vector<>::select_1_type sdb_sel(&classic);
    uint64_t a = sdb_sel(i_classic);
    uint64_t b;
    input >> b;//get first value for b
    while(i_classic <= nb_of_kmer && i_rc < nb_of_distinct_rev_comp){
        if(a < b){
            sdvB.set(a);
            i_classic++;
            a = sdb_sel(i_classic);
        } else if (b < a){
            sdvB.set(b);
            i_rc++;
            if(i_rc < nb_of_distinct_rev_comp) input >> b;//update value
        } else { // a = b
            cout << "Should not happen because we have removed the duplicate precedently" << endl;
        }
    }
    cout << "i_classic after while loop: " << i_classic << endl;
    cout << "i_rc after while loop: " << i_rc << endl;
    //at this point, either i_classic = classic.size() or i_rc = rc.size()
    if(i_classic <= nb_of_kmer){
        while(i_classic <= nb_of_kmer){
            sdvB.set(sdb_sel(i_classic));
            i_classic++;
        }
    } else { //i_rc < rc.size()
        while(i_rc < nb_of_distinct_rev_comp){
            sdvB.set(b);
            i_rc++;
            input >> b;
        }
    }
    auto e3 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed3 = e3 - b3;
    cout << "--> Time (s): " << elapsed3.count() << endl;

    cout << "i_classic at the end: " << i_classic << endl;
    cout << "i_rc at the end: " << i_rc << endl;
    input.close();
    sd_vector<> sdv(sdvB);
    return sdv;
}

/**
 * Returns a string representing the reverse complement of of 3-mer in forward.
 * @param seq
 * @return a string
 */
string decodeCaracterInForwardFormat(int seq){
    string res(3, ' ');
    //cout << "seq : " << seq << endl;
    for(int i(0); i < 3; i++){
        switch(seq & 0x3){ //compares the decimal value of the first two bits
            case 0: res[2-i] = 'T'; break;
            case 1: res[2-i] = 'G'; break;
            case 2: res[2-i] = 'C'; break;
            case 3: res[2-i] = 'A'; break;
        }
        seq >>= 2;
    }
    return res;
}

/**
 * Returns an int which represents the complement value of a nucleotid.
 * @param c
 * @return 0 or 1 or 2 or 3 according to c
 */
int encodeCaracterInReverseFormat(char c){
    int res = 0;
    switch(c){
        case 'A': res = 3; break;
        case 'C': res = 2; break;
        case 'G': res = 1; break;
        case 'T': res = 0; break;
    }
    return res;
}

//uncompress to classic version
/**
 * Uncompress a compressed version of a reverse complement.
 * @param compressedReverseComplement
 * @param sizeOfKmerInForward
 * @return a string which represents the forward version of a kmer.
 */
string uncompress(string compressedReverseComplement, int sizeOfKmerInForward){
    string res("");
    //cout << "compressedReverseComplement : " << compressedReverseComplement << endl;
    for(int i = compressedReverseComplement.size()-1; i >= 0; i--){
        string s = decodeCaracterInForwardFormat(compressedReverseComplement[i]-';');
        //cout << " s " << s << endl;
        res += s;
    }
    return res.substr(0, sizeOfKmerInForward);
}

/**
 * Compress a kmer.
 * @param compressedReverseComplement
 * @param sizeOfKmerInForward
 * @return a string which represents the compressed version of a kmer.
 */
string compress(string KmerForwardVersion){
    //int result_size = ceil(KmerForwardVersion.size()/3); //Returns the smallest integer that is greater than or equal to the float in paramater
    int nb_of_A_to_add = 0;
    while((KmerForwardVersion.size()+nb_of_A_to_add)%3 != 0) nb_of_A_to_add++;
    string s(nb_of_A_to_add,'A');
    //cout << "s size : " << s.size() << " s : " << s << endl;
    KmerForwardVersion += s;
    //cout << "KmerForwardVersion size : " << KmerForwardVersion.size() << " KmerForwardVersion : " << KmerForwardVersion << endl;
    int result_size = ceil(KmerForwardVersion.size()/3);
    string res(result_size, '*');
    int lastIndex = KmerForwardVersion.size()-1;
    for(int i = 0; i < result_size; i++){
        int referenceIndex = lastIndex-i*3;
        char firstCarac = 'A', secondCarac='A', thirdCarac='A';
        if(referenceIndex >= 0)   firstCarac  = KmerForwardVersion[referenceIndex];
        if(referenceIndex-1 >= 0) secondCarac = KmerForwardVersion[referenceIndex-1];
        if(referenceIndex-2 >= 0) thirdCarac  = KmerForwardVersion[referenceIndex-2];
        int trioValue = encodeCaracterInReverseFormat(firstCarac) +
        encodeCaracterInReverseFormat(secondCarac) * 4 + encodeCaracterInReverseFormat(thirdCarac) * 16;
        res[i] = ';'+trioValue;
    }
    return res;
}

/**
 * Returns an sd_vector representing a file's content. To do this, it stores in 
 * a external text file a compressed version of the reverse complements.
 * @param path - path of the file
 * @return an sd_vector representing file's content
 */
sd_vector<> fromFileToSdVector_TXTversionB64(string path){
    ifstream file(path, ios::in);  // Reading of the file which contains k-mers sequences
    if(!file) {
        cout << "Error while opening" << endl;
        return bit_vector{0};
    }

    //obtention of the k-mers' size
    string kmer("");
    file >> kmer;
    file.seekg(0, ios::beg);
    int K = kmer.size();
    
    //counting the number of k-mer in the file
    uint64_t nb_of_kmer = 0;
    string line("");
    while(getline(file, line)) nb_of_kmer++;
    cout << "nb_of_kmer : " << nb_of_kmer << endl;
    file.clear();
    file.seekg(0, ios::beg);
    sd_vector_builder sdv_builder_classic(pow(ALPHABET,K), nb_of_kmer);
    
    //first parsing to get the classic read
    while(file >> kmer){
        sdv_builder_classic.set(encode(kmer, K));//reverseComplementLexico(i_classic, K);
        file >> kmer;
    }
    file.clear();
    file.seekg(0, ios::beg);
    sd_vector<> classic(sdv_builder_classic);
    cout << "Length of classic : " << classic.size() << endl;
    cout << "size of classic (MB) : " << size_in_mega_bytes(classic) << endl;
    
    cout << "Completion of distinct_rev_comp.txt." << endl;
    auto b1 = std::chrono::high_resolution_clock::now();
    //second parsing to get the reverse complement
    ofstream output;
    output.open("distinct_rev_comp.txt");
    uint64_t nb_of_distinct_rev_comp = 0;
    while(file >> kmer){
        uint64_t i_classic = encode(kmer, K);
        uint64_t i_rc = reverseComplementLexico(i_classic, K);
        if(classic[i_rc]==0){ //if the reverse complement is not present in the classical reads then add it in the output
            //cout << "i : " << i_rc << endl;
            output << compress(kmer);//i_rc;
            output << '\n';
            nb_of_distinct_rev_comp++;
        }
        file >> kmer;
    }
    output.close();
    auto e1 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed1 = e1 - b1;
    cout << "--> Time (s): " << elapsed1.count() << endl;
    
    //sort the output
    cout << "distinct_rev_comp.txt completed. Now we sort." << endl;
    auto b2 = std::chrono::high_resolution_clock::now();
    string command = "sort distinct_rev_comp.txt -o distinct_rev_comp_sortedB64.txt";
    system(command.c_str());
    cout << "Sort finished. " << endl;
    auto e2 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed2 = e2 - b2;
    cout << "--> Time (s): " << elapsed2.count() << endl;
    cout << "Now we build the final sd_vector." << endl;
    
    //now we build the sd_vector
    auto b3 = std::chrono::high_resolution_clock::now();
    ifstream input;
    input.open("distinct_rev_comp_sortedB64.txt");  // Reading of the file which contains k-mers sequences
    if(!input) {
        cout << "Error while opening distinct_rev_comp.txt" << endl;
        return bit_vector{0};
    };
    sd_vector_builder sdvB(pow(ALPHABET,K), nb_of_kmer + nb_of_distinct_rev_comp);
    uint64_t i = 0, i_classic = 1, i_rc = 0, last_rc = 0;
    
    sd_vector<>::select_1_type sdb_sel(&classic);
    uint64_t a = sdb_sel(i_classic);
    uint64_t b;
    string b64;
    input >> b64;//get first value for b
    b = reverseComplementLexico(encode(uncompress(b64, K), K), K);
    while(i_classic <= nb_of_kmer && i_rc < nb_of_distinct_rev_comp){
        if(a < b){
            sdvB.set(a);
            i_classic++;
            a = sdb_sel(i_classic);
        } else if (b < a){
            sdvB.set(b);
            i_rc++;
            if(i_rc < nb_of_distinct_rev_comp){
                input >> b64;//get first value for b
                b = reverseComplementLexico(encode(uncompress(b64, K), K), K);
            }
        } else { // a = b
            cout << "Should not happen because we have removed the duplicate precedently" << endl;
        }
    }
    cout << "i_classic after while loop: " << i_classic << endl;
    cout << "i_rc after while loop: " << i_rc << endl;
    //at this point, either i_classic = classic.size() or i_rc = rc.size()
    if(i_classic <= nb_of_kmer){
        while(i_classic <= nb_of_kmer){
            sdvB.set(sdb_sel(i_classic));
            i_classic++;
        }
    } else { //i_rc < rc.size()
        while(i_rc < nb_of_distinct_rev_comp){
            sdvB.set(b);
            i_rc++;
            input >> b64;//get first value for b
            b = reverseComplementLexico(encode(uncompress(b64, K), K), K);
        }
    }
    auto e3 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed3 = e3 - b3;
    cout << "--> Time (s): " << elapsed3.count() << endl;
    
    cout << "i_classic at the end: " << i_classic << endl;
    cout << "i_rc at the end: " << i_rc << endl;
    input.close();
    sd_vector<> sdv(sdvB);
    return sdv;
}

