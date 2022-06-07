//
// Created by alexa on 18/06/2020.
//

/*
 * Return the merging version of 2 sd_vector
 * @param a an sd_vector of the same size as b
 * @param b an sd_vector
 * @param nb_of_1_in_A the number of 1-bit in the sd_vector a
 * @param nb_of_1_in_B the number of 1-bit in the sd_vector b
 * @return a sd_vector which is the merging version of the 2 sd_vector params
 */
sd_vector<> merge (const sd_vector<> &a, const sd_vector<> &b, int nb_of_1_in_A, int nb_of_1_in_B){
    int len = a.size();
    int nb_of_1 = nb_of_1_in_A + nb_of_1_in_B;
    int nb_of_1_in_merged = 0;
    sd_vector_builder s(len, nb_of_1);
    //filling builder
    for(int i = 0; i < len; i++) {
        if (a[i] == 1 && b[i] == 1)
            s.set(i);
        if (a[i] == 1 || b[i] == 1) {
            s.set(i);
            nb_of_1_in_merged++;
        }
    }
    //creation of the sd_vector
    sd_vector<> merger(s);
    return merger;
}

/**
 * Transform a number in base 10 in base 64.
 * @param valueInBase10
 * @return a string representing the number in base 64
 */
string convertToBase64(uint64_t valueInBase10){
    string res("");
    uint64_t biggest_divider = 1;
    while((uint64_t)(valueInBase10/biggest_divider) > 63) biggest_divider <<= 6;
    uint64_t remainder = valueInBase10;
    uint64_t quotient = valueInBase10;
    biggest_divider <<= 6; //*= 64 <->
    while(remainder != 0){
        biggest_divider >>= 6;
        quotient = remainder/biggest_divider;
        remainder = remainder % biggest_divider;
        res += (char)(';'+quotient);
    }
    if(biggest_divider > 1){
        string rest((int)(log(biggest_divider)/log(64)), ';');
        res += rest;
    }
    return res;
}

/**
 * Transform a number in base 10 in base 64.
 * @param valueInBase64
 * @return a number in base 10
 */
uint64_t convertFromBase64ToBase10(string valueInBase64){
    uint64_t res = 0;
    uint64_t pow = 1;
    for(int i = valueInBase64.size()-1; i >= 0; i--){
        res += (valueInBase64[i]-';') * pow;
        pow *= 64;
    }
    return res;
}

/* Calculate the reverse complement compressed version of a given k-mer
 * @param seq - the k-mer we want to know the reverse complement
 * @param len - the length of the total sequence of the given k-mer
 * @return a uint64_t which is the compressed version of the reverse complement
 */
uint64_t reverseComplement(string seq, uint64_t len){
    int sizeOfSeq = log(len) / log(ALPHABET);
    reverse(begin(seq), end(seq));  //Reverse of string
    uint64_t complem = encode(seq, sizeOfSeq);  //encoding the reverse
    uint64_t final;
    if(complem > len/2){    //Position for complement
        uint64_t position(len-complem);
        final = position-1;
    }else{
        uint64_t position(complem);
        final = len-position-1;
    }
    return final;
}

