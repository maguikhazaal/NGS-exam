# NGS-exam

On Gitbash 
First, we need to import the data from the shared data file /data-shared/vcf_examples/luscinia_vars.vcf.gz .  
We will store the file in a variable called _INPUT_ to make it easier to call it.

### Prepare the data
#### 1. Import data and define variable
```
INPUT=/data-shared/vcf_examples/luscinia_vars.vcf.gz
```
#### 2.Unzip data
The big data file is compressed in a vcf.gz format. We need to unzip it using the _zcat_ function. We can visualize it at any moment using _less -S_ function.  
```
<$INPUT zcat | less -S
```

####  3.Remove header:
The file is full of information that we do not necessarily need. We will progressively 'clean' the file to keep.
We will remove the header section, that is the one that starts with "##" , using the _grep -v_ function
```
<$INPUT zcat | grep -v "^##" | less -S
```

####  4.Filter out columns: '#CHROM', 'POS', 'QUAL':
For this task, we need the following variables only: chromosome name (#CHROM), position (POS), and PHRED score corresponding to quality (QUAL). We will filter out these variables that correspond to the columns 1,2, and 6, using the _cut_ function:
```
<$INPUT zcat | grep -v "^##" | cut -f1-2,6 | less -S
```

####  5.Get rid of the '#' before 'CHROM':
The '#' symbole before CHROM is not necessary, we will remove it as well, using the _tail_ function: 
```
<$INPUT zcat | grep -v "^##" | cut -f1-2,6 | tail -c +2 | less -S
```

####  6.Save data as tsv file
We will save the output file, that is, the "clean" data, in a _tab seperated variable_ format.
```
<$INPUT zcat | grep -v "^##" | cut -f1-2,4-6 | tail -c +2  > ~/projects/finalexam/data/luscinia_vars.tsv
```

#### 7. *_(OPTIONAL)_* Exclude the random , unmapped , and Ambiguous chromosomes:
We notice that for some chromosomes, they have the label 'random' (_random), unmapped (_Un ; _Unmapped), and 'ambiguous' (_Amb)._ We hypothesize that we can exclude these chromosomes from our analysis. We will remove it using the _grep -v_ function again, and store the new data file under the name 'luscinia_vars_select_chrom.tsv

```
<$INPUT zcat | grep -v "^##" | grep -v "random\>" | grep -v "Unmapped\>" | grep -v "Un\>" | grep -v "Amb\>" | cut -f1-2,6 | tail -c +2  > ~/projects/finalexam/data/luscinia_vars_select_chrom.tsv
```
