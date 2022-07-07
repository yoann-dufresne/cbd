#!/usr/bin/env nextflow
data=Channel.fromPath('./sorted/*')
rand = ['random','successive']
type = ['contain','successor']
percent = ['5','10','20','50','100']
number= ['1000','1000000','10000000']
process_number= Channel.of(1..10000)

data 
    .combine(number)
    .combine(rand)
    .combine(type)
    .set{data}


process test{
    publishDir "./result", mode: 'link' 
    input :
    tuple file(kmer),val(x),val(r),val(t) from data

    val pn from process_number
    output :
    file "result${pn}"
    script :
    """
    echo "${kmer}\nk=31\n${x}\n${r}\n${t}\n">result${pn}
    timerequest ${kmer} 31 ${x} ${r} ${t}  >> result${pn}
    """ 
}
process shuffle{
    input :
    file a from data2
    output :
    file "suffled_tmp" into shuf
    script:
    """
    shuf ${a} > shuffled_tmp
    """
}
data2=Channel.fromPath('./sorted/*')
type = ['contain','successor']
percent = ['5','10','20','50','100']
number= ['1000','1000000','10000000']
process_number= Channel.of(1..10000)

data2
    .combine(number)
    .combine(type)
    .combine(percent)
    .set{data2}


process percent {
    publishDir "./resultpercent", mode: 'link' 
    input :    
    tuple file(kmer),val(x),val(t),val(percent) from data2

    file shuf from shuf

    val pn from process_number
    output :
    file "result${pn}"

    script:
    """
    echo "${kmer}\nk=31\n${x}\nrandom\n${t}\n${percent}">result${pn}
    timerequest ${kmer} 31 ${x} random ${t} ${percent} ${shuf} >> result${pn}

    """


}



