#Make script executable using 
chmod +x unix_final_exam_workflow.sh

#1. Optional, create directory. If needed, remove the '#' symbole before the next line
#-------------------------------------
#mkdir -p projects/unix_final_exam/data

#2. Import data and assign variable "INPUT"
#------------------------------------------
INPUT=/data-shared/vcf_examples/luscinia_vars.vcf.gz
echo "Data stored in variable INPUT"

#3. Unzip data file using zcat, remove header using grep -v, filter out chosen columns (CHROM, POS, QUAL) corresponding to columns 1,2 and 6, using cut, and remove '#' symbole before CHROM using tail.
#Save resulting output as a tab seperated variable file. 
#------------------------------------------------------------
<$INPUT zcat | grep -v "^##" | cut -f1-2,6 | tail -c +2  > ~/projects/unix_final_exam/data/luscinia_vars.tsv
echo "Output : filtered columns inc all chrom > luscinia_vars.tsv"

#4. Filter out chromosomes labeled 'Random', 'Unmapped', or 'Ambiguous' from file. Save output as luscinia_vars_select_chrom.tsv
#----------------------------------------------------------------------------------------------------------------------------------
<$INPUT zcat | grep -v "^##" | grep -v "random\>" | grep -v "Unmapped\>" | grep -v "Un\>" | grep -v "Amb\>" | cut -f1-2,6 | tail -c +2  > ~/projects/unix_final_exam/data/luscinia_vars_select_chrom.tsv
echo "Output : filtered columns exc random, unmapped, and ambiguous chrom > luscinia_vars_select_chrom.tsv"


#5. Echo "Done!" and move to R processing  
#-----------------------------------------
echo "Done!"
echo "Move to R"

Rscript plot_final_exam/unix_final_exam_Rscript.R

