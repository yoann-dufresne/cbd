#!/usr/bin/env nextflow
data=Channel.fromPath('./sorted/*')
rand = ['0','1']
type = ['contain','successor']
library = ['BM','SDSL']

process testcontainsBM{
    publishDir "./result", mode: 'link' 
    input :
    file kmer from data
    each r from rand
    each t from type
    each lib from library
    output :
    file "result_${r}_${t}_${lib}"
    script :
    """
    main ${kmer} 31 10 ${r} ${t} ${lib} > result_${r}_${t}_${lib}
    """ 
}

