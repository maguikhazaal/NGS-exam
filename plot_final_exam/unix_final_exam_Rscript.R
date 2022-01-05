setwd("/home/user01/projects/unix_final_exam/plot_final_exam/")

# 1. Prepare library
#================================
library (tidyverse)
library (dplyr)
library (ggplot2)
# 2. Prepare data
#=================================

# 2.1) Load input
#----------------

read_tsv('~/projects/unix_final_exam/data/luscinia_vars_select_chrom.tsv') -> data_select_chrom
#View (data_select_chrom)


# 2.2) Optional, Remove white space in POS and QUAL (for aesthetics)
#-----------------------------------------

data_select_chrom$POS <- trimws(data_select_chrom$POS, which = c("both"))
data_select_chrom$QUAL <- trimws(data_select_chrom$QUAL, which = c("both"))

# 2.3) Mutate CHROM from integer to factor
#-------------------------------------------

read_tsv('~/projects/finalexam/data/luscinia_vars_norandom.tsv') %>%
  mutate(CHROM = as.factor(CHROM)) ->
  data_select_chrom

# 3. Plots
#===============================================
#3.1 Distribution of PHRED Quality over whole Genome
#---------------------------------------------------
p1 = data_select_chrom %>%
  ggplot(aes(QUAL)) +
  geom_histogram(binwidth=1) +
  ylab("Count of variants") +
  xlab("PHRED quality") +
  ggtitle("Distribution of PHRED score quality on whole genome wrt to Count of Variants")

plot (p1)

#We notice that that some observations have a PHRED quality score of 999, which probably corresponds to an error of some sort. We will consider these values NA and exclude them from our analysis. Each plot will be saved a a PDF file.

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


# For added information, we can group the chromosomes by colors

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


#This is not very informative. We will represent the phred score distribution based on the position in the genome


#data_select_chrom %>% 
#  filter(QUAL != 999) %>% 
#  ggplot(aes(POS, QUAL, colour = CHROM)) + 
#  geom_point(size = 0.05) +
#  theme(legend.key.height= unit(0.30, 'cm'), legend.key.width= unit(1, 'cm'), legend.key.size = unit(1, 'cm')) +
#  xlab('POSITION') + 
#  ylab('QUALITY') + 
#  ggtitle("Distribution of PHRED quality over whole genome and by chrom")

#This representation is too noisy, we will use a boxplot instead

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

##To make it even clearer, we will define background colors depending on the PHRED score range, as follows: Zone1 : score 0-20 (yellow) Zone2 : score 20 - 60 (green) Zone3 : score 60 - 400 (red)
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
  

  #3.2 Distribution of Quality on Chromosome
  #-------------------------------------------
  
  ## Alternatively, we can represent the PHRED quality distribution on each chromosome 
  ## Using a dot plot 
#  data_select_chrom %>% 
#    filter(QUAL != 999) %>% 
#    ggplot() + 
#    geom_point(aes(POS, QUAL, colour = CHROM), size = 0.05) +
#    theme(legend.key.height= unit(0.5, 'cm'), legend.key.width= unit(1, 'cm'), legend.key.size = unit(1, 'cm')) +
#    xlab('POSITION') + 
#    ylab('QUALITY') + 
#    ggtitle("Distribution of PHRED quality by chrom")
 
   ## This is rather noisy, better representation would be a box plot
  
  ##As a box plot (final)
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





savehistory()
