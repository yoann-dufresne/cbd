#!/usr/bin/env nextflow
data=Channel.fromPath('./sorted/*')
rand = ['random','successive']
type = ['contain','successor']
library = ['BM','SDSL']
number= ['1','10','100','1000','1000000','1000000000']
process_number= Channel.of(1..10000)
data 
    .combine(number)
    .combine(rand)
    .combine(type)
    .combine(library)
    .set{data}
process test{
    publishDir "./result", mode: 'link' 
    input :
    tuple file(kmer),val(x),val(r),val(t),val(lib) from data

    val pn from process_number
    output :
    file "result${pn}"
    script :
    """
    echo "${kmer}\nk=31\n${x}\n${r}\n${t}\n${lib}">result${pn}
    timerequest ${kmer} 31 ${x} ${r} ${t} ${lib} >> result${pn}
    """ 
}

