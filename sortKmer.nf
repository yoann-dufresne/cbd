#!/usr/bin/env nextflow

data=Channel.fromPath('./data/*.fasta')


process countKmer {         
    input :
    path dataentry from data
    
    output :
    tuple file("${dataentry.baseName}.kmc_pre"), file("${dataentry.baseName}.kmc_suf") into bin

    script:
    """
    mkdir bloup
    kmc -fm -ci1 ${dataentry} ${dataentry.baseName} .
    """
}

process convert {
    input :                                                                      
    tuple file(pre), file(suf) from bin
    output :
    file("dataout/${pre.baseName}") into txt

    script :
    """
    mkdir dataout
    kmc_tools transform  ${pre.baseName} dump dataout/${pre.baseName}
    """
}

//sort les kmer après qu'ils ait été mis en forme 
process sortKmer { 
    publishDir "./sorted/", mode: 'link' 
                                                                                                                   
    input:
    file convert from txt
    output:
    file "sorted_${convert.baseName}"
    
    """
    sort ${convert} > sorted_${convert.baseName}
    """


}
