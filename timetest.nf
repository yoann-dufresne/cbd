#!/usr/bin/env nextflow
kmer=Channel.fromPath('./sorted/*')
process serialize{
    storeDir'./serial'
    input:
    file(kmer) from kmer
    output:
    tuple file("serial${kmer}"),file(kmer) into serial
    script:
    """
    serializer ${kmer} serial${kmer}
    """
}

type = ['contains','successor']
iteration=Channel.of(1..10)
serial
    .combine(type) 
    .combine(iteration)
    .into{data;data2;data3}

rand = ['random','sequence']
number= ['1000']
data 
    .combine(rand)
    .combine(number)
    .set{data}


process_number= Channel.of(1..10000)
process test{
    publishDir "./result", mode: 'link' 
    input :
    tuple file(serial),file(kmer),val(t),val(i),val(r),val(x) from data
    val pn from process_number
    output :
    file "result${pn}"

    memory 1.GB
    cpus 1
    script :
    """
    echo "${kmer}\nk=31\n${x}\n${r}\n${t}\n${i}\n">result${pn}
    timerequest ${serial} 31 ${x} ${r} ${t}  >> result${pn}
    """ 
}





percent = ['0','5','10','15','20','25','30','35','40','45','50','55','60','65','70','75','80','85','90','95','100']
number= ['100000']
data2
    .combine(percent)
    .combine(number)
    .set{data2}

process_number= Channel.of(1..10000)
process percent {
    publishDir "./resultpercent", mode: 'link' 
    input :    
    tuple file(serial),file(kmer),val(t),val(i),val(percent),val(x) from data2
    val pn from process_number
    output :
    file "result${pn}"

    memory 1.GB
    cpus 1
    script:
    """
    shuf ${kmer} -n ${x} > tmpshuf
    echo "${kmer}\nk=31\n${x}\nrandom\n${t}\n${percent}">result${pn}
    timerequest ${serial} 31 ${x} random ${t} ${percent} tmpshuf >> result${pn}
    echo "${i}">>result${pn}
    """


}
percent = ['0','5','10','15','20','25','30','35','40','45','50','55','60','65','70','75','80','85','90','95','100']
number= ['1000']
fasta=Channel.fromPath("./data/*")
shufsequence=Channel.fromPath('./sequencegen.py')
data3
    .combine(percent)
    .combine(number)
    .combine(fasta)
    .combine(shufsequence)
    .set{data3}

process_number= Channel.of(1..10000)
process percentsequence{
    publishDir "./resultpercentseq", mode: 'link'
    input:
    tuple file(serial),file(kmer),val(t),val(i),val(percent),val(x),file(fasta),file(shufsequence) from data3
    val pn from process_number

    output:
    file "result${pn}"

    memory 1.GB
    cpus 1
    script:
    """
    python3 ${shufsequence} ${fasta} >shuf
    echo "${kmer}\nk=31\n${x}\nsequence\n${t}\n${percent}">result${pn}
    timerequest ${serial} 31 ${x} sequence ${t} ${percent} shuf >> result${pn}
    echo "${i}">>result${pn}
    """
}














