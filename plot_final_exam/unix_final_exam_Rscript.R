setwd("/home/user01/projects/unix_final_exam/plot_final_exam")

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
data_select_chrom %>% 
  ggplot(aes(QUAL)) +
  geom_histogram(binwidth=1) +
  ylab("Count of variants") +
  xlab("PHRED quality") +
  ggtitle("Distribution of PHRED score quality on whole genome wrt to Count of Variants")

#We notice that some observations have a PHRED score of 999. We will consider these values to be NA and exclude them.

data_select_chrom %>% 
  filter(QUAL < 999) %>% 
  ggplot(aes(QUAL)) +
  geom_histogram(binwidth=1) +
  ylab("Count of variants") +
  xlab("PHRED quality") +
  ggtitle("Distribution of PHRED score quality on whole genome wrt Count of Variants")

# For added information, we can group the chromosomes by colors

data_select_chrom %>% 
  filter(QUAL != 999) %>% 
  ggplot(aes(QUAL, colour = CHROM, group = CHROM)) +
  geom_histogram(binwidth=1) +
  ylab("Count of variants") +
  xlab("PHRED quality") +
  ggtitle("Distribution of PHRED score quality on whole genome wrt Count of Variants") +
  theme(legend.key.height= unit(0.25, 'cm'), legend.key.width= unit(0.25, 'cm'), legend.key.size = unit(0.5, 'cm'))

#This is not very informative. We will represent the phred score distribution based on the position in the genome
data_select_chrom %>% 
  filter(QUAL != 999) %>% 
  ggplot(aes(POS, QUAL, colour = CHROM)) + 
  geom_point(size = 0.05) +
  theme(legend.key.height= unit(0.30, 'cm'), legend.key.width= unit(1, 'cm'), legend.key.size = unit(1, 'cm')) +
  xlab('POSITION') + 
  ylab('QUALITY') + 
  ggtitle("Distribution of PHRED quality over whole genome and by chrom")

#Again, this is not the best representation. We will try to represent in a boxplot

data_select_chrom %>% 
  filter(QUAL != 999) %>% 
  ggplot(aes(POS, QUAL, fill = CHROM)) + 
  geom_boxplot() +
  theme(legend.key.height= unit(0.25, 'cm'), legend.key.width= unit(0.25, 'cm'), legend.key.size = unit(1, 'cm')) +
  xlab('POSITION') + 
  ylab('QUALITY (PHRED Score)') + 
  scale_y_log10() +
  scale_x_log10() +
  ggtitle("Distribution of PHRED quality over Whole Genome and by Chromosome")

##Define zones
data.frame(
  ymin = c(0, 20, 60),
  ymax = c(20, 60, 400),
  Colour = c("yellow", "lawngreen", "violetred2")) ->
  quals


  data_select_chrom %>% 
  filter(QUAL != 999) %>% 
  ggplot() +
    geom_rect(aes(ymin = ymin, ymax = ymax, fill = Colour),
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
    ggtitle("Distribution of PHRED quality over Whole Genome and by Chromosome")

  #3.2 Distribution of Quality on Chromosome
  #-------------------------------------------
  
  ## As dot plot (final)
  data_select_chrom %>% 
    filter(QUAL != 999) %>% 
    ggplot() + 
    geom_point(aes(POS, QUAL, colour = CHROM), size = 0.05) +
    theme(legend.key.height= unit(0.5, 'cm'), legend.key.width= unit(1, 'cm'), legend.key.size = unit(1, 'cm')) +
    xlab('POSITION') + 
    ylab('QUALITY') + 
    ggtitle("Distribution of PHRED quality by chrom")
  ## This is rather noisy, better representation would be a box plot
  
  ##As a box plot (final)
  data_select_chrom %>% 
    filter(QUAL != 999) %>%
    ggplot() +
    geom_rect(aes(ymin = ymin, ymax = ymax, fill = Colour),
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
  
  
  savehistory()
