#!/usr/bin/env nextflow
data=Channel.fromPath('./sorted/*')
rand = ['random','sequence']
type = ['contain','successor']
percent = ['5','10','20','50','100']
number= ['1000']
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

    cpus 1
    memory 13.GB
    script :
    """
    echo "${kmer}\nk=31\n${x}\n${r}\n${t}\n">result${pn}
    timerequest ${kmer} 31 ${x} ${r} ${t}  >> result${pn}
    """ 
}
data4=Channel.fromPath('./sorted/*')
process shuffle{
    input :
    file a from data4
    output :
    file "shuffled_tmp" into shuf
    script:
    """
    shuf ${a} > shuffled_tmp
    """
}
data2=Channel.fromPath('./sorted/*')
type = ['contain','successor']
percent = ['5','10','20','50','100']
number= ['1000']
process_number= Channel.of(1..10000)

data2
    .combine(number)
    .combine(type)
    .combine(percent)
    .combine(shuf)
    .set{data2}


process percent {
    publishDir "./resultpercentseq", mode: 'link' 
    input :    
    tuple file(kmer),val(x),val(t),val(percent),file(shuf) from data2
    val pn from process_number
    output :
    file "result${pn}"

    script:
    """
    echo "${kmer}\nk=31\n${x}\nrandom\n${t}\n${percent}%">result${pn}
    timerequest ${kmer} 31 ${x} random ${t} ${percent} ${shuf} >> result${pn}

    """


}
data3=Channel.fromPath('./sorted/*')
type = ['contain','successor']
percent = ['5','10','20','50','100']
number= ['10']
process_number= Channel.of(1..10000)
sequences=Channel.fromPath('./sequences/*')
data3
    .combine(number)
    .combine(type)
    .combine(percent)
    .combine(sequences)
    .set{data3}


process percentsequence{
    publishDir "./resultpercent", mode: 'link'
    input:
    tuple file(kmer),val(x),val(t),val(percent),val(s) from data3
    val pn from process_number

    output:
    file "result${pn}"

    script:
    """
    echo "${kmer}\nk=31\n${x}\nrandom\n${t}\n${percent}%">result${pn}
    timerequest ${kmer} 31 ${x} sequence ${t} ${percent} ${s} >> result${pn}
    """



}



