#!/usr/bin/env nextflow
data=Channel.fromPath('./sorted/*')
rand = ['random','sequence']
type = ['contain','successor']
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

    memory 15.GB
    cpus 1
    script :
    """
    echo "${kmer}\nk=31\n${x}\n${r}\n${t}\n">result${pn}
    timerequest ${kmer} 31 ${x} ${r} ${t}  >> result${pn}
    """ 
}

data2=Channel.fromPath('./sorted/*')
type = ['contain','successor']
percent = ['5','10','15','20','25','30','35','40','45','50','55','60','65','70','75','80','85','90','95','100']
number= ['100000']
iteration=channel.of(1..10)
process_number= Channel.of(1..10000)

data2
    .combine(number)
    .combine(type)
    .combine(percent)
    .combine(iteration)
    .set{data2}


process percent {
    publishDir "./resultpercent", mode: 'link' 
    input :    
    tuple file(kmer),val(x),val(t),val(percent),val(i) from data2
    val pn from process_number
    output :
    file "result${pn}"

    memory 15.GB
    cpus 1
    script:
    """
    shuf ${kmer} -n ${x} > tmpshuf
    echo "${kmer}\nk=31\n${x}\nrandom\n${t}\n${percent}">result${pn}
    timerequest ${kmer} 31 ${x} random ${t} ${percent} tmpshuf >> result${pn}
    echo "${i}">>result${pn}
    """


}
data3=Channel.fromPath('./sorted/*')
type = ['contain','successor']
percent = ['5','10','15','20','25','30','35','40','45','50','55','60','65','70','75','80','85','90','95','100']
number= ['1000']
process_number= Channel.of(1..10000)
iteration=Channel.of(1..10)
shufsequence=Channel.fromPath('./sequencegen.py')
fasta=Channel.fromPath("./data")
data3
    .combine(number)
    .combine(type)
    .combine(percent)
    .combine(iteration)
    .combine(shufsequence)
    .combine(fasta)
    .set{data3}

process percentsequence{
    publishDir "./resultpercentseq", mode: 'link'
    input:
    tuple file(kmer),val(x),val(t),val(percent),val(i),file(shufsequence),file(fasta) from data3
    val pn from process_number

    output:
    file "result${pn}"

    memory 15.GB
    cpus 1
    script:
    """
    python3 ${shufsequence} ${fasta} >shuf
    echo "${kmer}\nk=31\n${x}\nsequence\n${t}\n${percent}">result${pn}
    timerequest ${kmer} 31 ${x} sequence ${t} ${percent} shuf >> result${pn}
    echo "${i}">>result${pn}
    """



}



