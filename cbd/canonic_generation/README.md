# Scripts for helping units tests
## involveGenerate_counts.py
Python script which is made to generate canonical k-mers only.<br>
3 options :<br>
**-k** : the size of each k-mers<br>
**-n** : the length of the nucleotides sequence.<br>
**-f** : encoding choice, can be ACGT or ACTG.<br>
Examples : 
```
#I want to print 100 4-mers canonical for ACTG :
python3 involveGenerate_counts.py -k 4 -n 100 -f ACTG
```
```
#I want to write the information above in a file : 
python3 involveGenerate_counts.py -k 4 -n 100 -f ACTG > unsortACTG.txt
```
## sortACTG.py
Python script to help sorting elements of a file in ACTG encoding.<br>
We use this special sort because the classic shell sort is made to sort string lexicographicaly.<br>
ACTG format can't be sort lexicographically due to restrictions on **fromFileToSdVector**.<br>
If we want to sort an ACGT encoding file, just need shell sort : 
``
sort unsortACGT.txt > sortACGT.txt
``
<br>Example :
```
#I want to sort a file with ACTG encoding :
python 3 sortACTG.py unsortACTG.py > sortACTG.py 
```
