#!/usr/bin/env nextflow
data=Channel.fromPath('./sorted/*')
rand = ['random','successive']
type = ['contain','successor']
library = ['SDSL']
number= ['1000','1000000','1000000000']
iteration=Channel.of(1..10)
process_number= Channel.of(1..10000)
data 
    .combine(number)
    .combine(rand)
    .combine(type)
    .combine(library)
    .combine(iteration)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
    .set{data}
process test{
    publishDir "./result", mode: 'link' 
    input :
    tuple file(kmer),val(x),val(r),val(t),val(lib),val(i) from data

    val pn from process_number
    output :
    file "result${pn}"

    cpus 1
    memory 13.GB
    script :
    """
    echo "${kmer}\nk=31\n${x}\n${r}\n${t}\n${lib}\n iteration number ${i}">result${pn}
    timerequest ${kmer} 31 ${x} ${r} ${t} ${lib} >> result${pn}
    """ 
}

