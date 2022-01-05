# MB170C47: UNIX and work with genomic data
The present repository contains one solution to the final task in the MB170C47 course: "Unix and Work with Genomic Data" from the Faculty of Science, Charles University, Prague. 
This project is based on the sample vcf file /luscinia_vars.vcf.gz available on the shared data prepared for the course. The data file contains genomic information about the luscinia species that can be used in various analyses. In this part, we will be analyzing the distribution of the PHRED qualities on the whole genome and by chromosome. 
The project will be carried out in two parts: processing the vcf file in Unix, and a graphical representation in R. 

# PROCESSING DATA IN UNIX

### I. Prepare working space
#### 1. Create directory 
We will create a directory 'unix_final_exam' and work in the data file 
```
mkdir -p projects/unix_final_exam/data
```

#### 2. Note and save working code 

```
cd projects/unix_final_exam
nano unix_final_exam_workflow.sh
```

#### 3. Make script executable
```
chmod +x unix_final_exam_workflow.sh
```

### II. Prepare the data
#### 1. Import data and define variable
First, we need to import the data from the shared data file /data-shared/vcf_examples/luscinia_vars.vcf.gz .  
We will store the file in a variable called _INPUT_ to make it easier to call it.
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
<$INPUT zcat | grep -v "^##" | cut -f1-2,6 | tail -c +2  > ~/projects/unix_final_exam/data/luscinia_vars.tsv
```

#### 7. *_(OPTIONAL)_* Exclude the random , unmapped , and Ambiguous chromosomes:
We notice that for some chromosomes, they have the label 'random' (_random), unmapped (_Un ; _Unmapped), and 'ambiguous' (_Amb)._ We hypothesize that we can exclude these chromosomes from our analysis. We will remove it using the _grep -v_ function again, and store the new data file under the name 'luscinia_vars_select_chrom.tsv

```
<$INPUT zcat | grep -v "^##" | grep -v "random\>" | grep -v "Unmapped\>" | grep -v "Un\>" | grep -v "Amb\>" | cut -f1-2,6 | tail -c +2  > ~/projects/unix_final_exam/data/luscinia_vars_select_chrom.tsv
```

---
# PROCESSING DATA IN R 
### I. (Optional) Set working directory
setwd("/home/user01/projects/unix_final_exam/plot_final_exam")

### II. Prepare the data

#### 1. Prepare library 
```
library (tidyverse)
library (dplyr)
library (ggplot)
```

### 2. Prepare data
#### 2.1) Load input
```
read_tsv('~/projects/unix_final_exam/data/luscinia_vars_select_chrom.tsv') -> data_select_chrom
View (data_select_chrom)
```

#### 2.2) Optional, Remove white space in POS and QUAL (for aesthetics)
```
data_select_chrom$POS <- trimws(data_select_chrom$POS, which = c("both"))
data_select_chrom$QUAL <- trimws(data_select_chrom$QUAL, which = c("both"))
```

#### 2.3) Mutate CHROM from integer to factor
```
read_tsv('~/projects/finalexam/data/luscinia_vars_norandom.tsv') %>%
  mutate(CHROM = as.factor(CHROM)) ->
  data_select_chrom
```

### 3. Plots

#### 3.1 Distribution of PHRED Quality over whole Genome
```
p1 = data_select_chrom %>%
  ggplot(aes(QUAL)) +
  geom_histogram(binwidth=1) +
  ylab("Count of variants") +
  xlab("PHRED quality") +
  ggtitle("Distribution of PHRED score quality on whole genome wrt to Count of Variants")

plot (p1)
```

We notice that that some observations have a PHRED quality score of 999, which probably corresponds to an error of some sort. We will consider these values NA and exclude them from our analysis. 
Each plot will be saved a a PDF file. 

```
p2 = data_select_chrom %>% 
  filter(QUAL != 999) %>% 
  ggplot(aes(QUAL)) +
  geom_histogram(binwidth=1) +
  ylab("Count of variants") +
  xlab("PHRED quality") +
  ggtitle("Distribution of PHRED score quality on whole genome wrt Count of Variants")

pdf("results/phred-genome-cov.pdf",w=10,h=8)
plot(p2)
dev.off()
```
![image](https://user-images.githubusercontent.com/83076900/148205922-673610a3-1e10-43f0-808e-0f9bc64828df.png)


For added information, we can group the chromosomes by colors

```
p3 = data_select_chrom %>% 
  filter(QUAL != 999) %>% 
  ggplot(aes(QUAL, colour = CHROM, group = CHROM)) +
  geom_histogram(binwidth=1) +
  ylab("Count of variants") +
  xlab("PHRED quality") +
  ggtitle("Distribution of PHRED score quality on whole genome wrt Count of Variants by chrom") +
  theme(legend.key.height= unit(0.25, 'cm'), legend.key.width= unit(0.25, 'cm'), legend.key.size = unit(0.5, 'cm'))

pdf("results/phred-genome-cov-by-chrom.pdf",w=10,h=8)
plot(p3)
dev.off()

```

![image](https://user-images.githubusercontent.com/83076900/148206051-c4b8457e-b1e0-49bc-bd8c-9f438cea593b.png)

This representation however is not very informative. We will represent the PHRED quality distribution on the whole genome and by chromosome based on the bp position. 

```
p4 = data_select_chrom %>% 
  filter(QUAL != 999) %>% 
  ggplot(aes(POS, QUAL, colour = CHROM)) + 
  geom_boxplot() +
  theme(legend.key.height= unit(0.25, 'cm'), legend.key.width= unit(0.25, 'cm'), legend.key.size = unit(1, 'cm')) +
  xlab('POSITION') + 
  ylab('QUALITY (PHRED Score)') + 
  scale_y_log10() +
  scale_x_log10() +
  ggtitle("Distribution of PHRED quality over Whole Genome and by Chromosome wrt Position")

pdf("results/phred-genome-pos-by-chrom.pdf",w=10,h=8)
plot(p4)
dev.off()


 ```
![image](https://user-images.githubusercontent.com/83076900/148206280-f64b5174-06c7-41a1-a509-41b9fd469e18.png)


To make it even clearer, we will define background colors depending on the PHRED score range, as follows: 
Zone1 : score 0-20 (yellow)
Zone2 : score 20 - 60 (green)
Zone3 : score 60 - 400 (red) 

```
data.frame(
  ymin = c(0, 20, 60),
  ymax = c(20, 60, 400),
  colour = c("yellow", "lawngreen", "violetred2")) ->
  quals


  p5 = data_select_chrom %>% 
  filter(QUAL != 999) %>% 
  ggplot() +
    geom_rect(aes(ymin = ymin, ymax = ymax, fill = colour),
              xmin = -Inf,
              xmax = Inf,
              alpha=0.3,
              data = quals,
              show.legend = FALSE) +
    scale_fill_identity() +
    geom_boxplot(aes(POS, QUAL, colour = CHROM), outlier.colour = NA) +
    theme(legend.key.height= unit(0.25, 'cm'), legend.key.width= unit(0.25, 'cm'), legend.key.size = unit(1, 'cm')) +
    xlab('POSITION') + 
    ylab('QUALITY (PHRED Score)') + 
    scale_y_log10() +
    scale_x_log10() +
    ggtitle("Distribution of PHRED quality over Whole Genome and by Chromosome wrt Position")
  
  pdf("results/phred-genome-pos-by-chrom-wzone.pdf",w=10,h=8)
  plot(p5)
  dev.off()


 ```
![image](https://user-images.githubusercontent.com/83076900/148206541-fc709918-c7bf-4298-912f-823e0f01dadc.png)



Alternatively, we can represent the PHRED quality distribution on each chromosome using a boxplot

```
  p6 = data_select_chrom %>% 
    filter(QUAL != 999) %>%
    ggplot() +
    geom_rect(aes(ymin = ymin, ymax = ymax, fill = colour),
              xmin = -Inf,
              xmax = Inf,
              alpha=0.5,
              data = quals,
              show.legend = FALSE) +
    scale_fill_identity() +
    geom_boxplot(aes(CHROM, QUAL, fill = "grey78")) +
    geom_smooth(aes(CHROM, QUAL, group = 1), colour = "blue") +
    scale_y_log10() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.key.size = unit(0.25, 'cm')) +
    xlab("Chromosome")+
    ylab("PHRED Quality") + 
    ggtitle("Distribution of PHRED quality by chromosome")
  
pdf("results/phred-by-chrom.pdf",w=10,h=8)
plot(p6)
dev.off()

```

![image](https://user-images.githubusercontent.com/83076900/148206663-85ed34f1-1cb1-4002-acef-fa03e81e50a1.png)
