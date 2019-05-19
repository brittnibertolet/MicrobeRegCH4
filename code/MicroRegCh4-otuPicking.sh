#################################################
# Bash script for OTU picking using QIIME 
# Raw sequence reads are available online at DOE JGI Genome Portal under Project ID 1041357. 
#################################################

#### Demultiplex reads and quality filter
split_libraries_fastq.py -o methanogenSurvey/ -i forward_reads.fastq.gz -b barcodes.fastq.gz -m map.tsv

#### Pick OTUs using open reference picking 
pick_open_reference_otus.py -o qiime_otus/ -i methanogenSurvey/seqs.fna -p ../uc_fast_params.txt

#### Align representative sequences with default pynast method
align_seqs.py -i rep_set.fasta -t core_set_template.fasta -o pynast_aligned/

### Assign taxonomy 
assign_taxonomy.py -i rep_set_aligned.fasta -m rdp

#### Check for chimeras 
# Identify chimeras from representative sequences
identify_chimeric_seqs.py -m ChimeraSlayer -i rep_set_aligned.fasta -a reference_set_aligned.fasta -o chimeric_seqs.txt
# Remove chimeras from representative sequences
filter_fasta.py -f rep_set_aligned.fasta -o non_chimeric_rep_set_aligned.fasta -s chimeric_seqs.txt -n
# Remove chimeras from OTU table
make_otu_table.py -i otu_map.txt -o otu_table.biom -e chimeric_seqs.txt -t taxonomy.txt

