#!/usr/bin/env nextflow

data=Channel.fromPath('/home/oceane/dev/data/*.fasta')


process countKmer {         
    input :
    path dataentry from data
    
    output :
    tuple file("${dataentry.baseName}.kmc_pre"), file("${dataentry.baseName}.kmc_suf") into bin
    
    memory 5.GB

    script:
    """
    mkdir bloup
    kmc -k31 -fm -ci1 ${dataentry} ${dataentry.baseName} .
    """
}

process convert {
    input :                                                                      
    tuple file(pre), file(suf) from bin
    output :
    file("dataout/${pre.baseName}") into txt

    memory 5.GB

    script :
    """
    mkdir dataout
    kmc_tools transform  ${pre.baseName} dump dataout/${pre.baseName}
    """
}

//sort les kmer après qu'ils ait été mis en forme 
process sortKmer { 
    publishDir "/home/oceane/dev/sorted/", mode: 'link' 
                                                                                                                   
    input:
    file convert from txt
    output:
    file "sorted_${convert.baseName}"

    memory 5.GB
    cpus 4
    
    """
    sort --parallel=4 ${convert} > sorted_${convert.baseName}
    """


}
