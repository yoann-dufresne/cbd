#!/usr/bin/env nextflow

data=Channel.fromPath('./data/*.fasta')

data
    .view()
    .into{d1;d2}

process countKmer {
    input :
    path dataentry from d1
    /*
    output :
    file("${dataentry}.tmp") into bin


    script:*/
    """
    kmc -fm -ci1 /home/oceane/dev/data/unknown.fasta /home/oceane/dev/data/bloup .
    """
}
// bin
//     .view()
//     .into{bin2}
// process convert {
//     input :
//     path dataentry from bin2
//     script :
//     """
//     kmc_tools transform  ${dataentry}.tmp dump test
//     """
// }
/*
//sort les kmer après qu'ils ait été mis en forme 
process sortKmer {
    input:
    path unsorted from counted
    output:
    
    """
    ./sort ${unsorted} > sorted
    """


}

process timeCBD {


}

process lineSplitter {
    input:


}*/