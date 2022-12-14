# Construct-ancestral-genome

Construct the ancestral genome from the given population(s).

https://github.com/Shuhua-Group/Construct-ancestral-genome

https://github.com/zhangxx123456/Construct-ancestral-genome




## Ussage:

python3 construct_ancestral_genome.py \\ <br>Outputfile  \\ <br>     Proportion_cutoff_of_haplotypes   \\ <br>    Recombination_distance  \\ <br>   Inputfiles  



## Parameters:

***Outputfile***: The output file for recording results of a constructed ancestral haplotype.

***Proportion_cutoff_of_haplotypes***: The lowest proportion of inferred ancestral haplotypes in all of the input files. Range from 0.0-1.0. For example, 0.2 means more than 20% of the haplotypes from input files in the physical position are available to provide gene pool for further construction.

***Recombination_distance***:  If the selected haplotype was extended to the given length, it will occur a recombination event. The unit of parameter is basepair(bp).

***Inputfiles***: One file or multiple files are allowed. All of the input files must have the completely consistent bi-allelic variants. '0' means the reference allele was selected to construct ancestral genomes. '1' means the alternative allele was selected to construct ancestral genomes. '.' means the allele was selected to construct ancestral genomes.



## For example:
python3 construct_ancestral_genome.py  \\ <br>   abc.chr1.txt.gz   \\ <br>  0.2   \\ <br>  6000  \\ <br>   Test.dataset/A.chr1.txt.gz   Test.dataset/B.chr1.txt.gz   Test.dataset/C.chr1.txt.gz


