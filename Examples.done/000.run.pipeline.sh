# This pipeline performs local ancestry inference for target populations using ChromoPainterV2. 
# It processes genomic data, infers local ancestry, and assembles ancestral haplotypes. 
# 
# Author: Xiaoxi Zhang
# Purpose: To process genomic data for local ancestry inference and assemble ancestry-specific haplotypes.
# To help users run the pipeline more efficiently, we have removed the data preprocessing section and the subsequent analysis for generating samples. If users require these steps, please refer to our paper.
# 
# Pipeline Steps:
# 1. Local ancestry inference for the target populations.
# 2. Determine which ancestries the recipient haplotypes copies from at each allele.
# 3. Data transpose and adding physical location information.
# 4. Extract the specified ancestral alleles from the selected samples.
# 5. Assemble ancestral haplotypes using extracted ancestral gene pools.
# 6. Generate ancestral individuals using assembled ancestral haplotypes.
# 7. Synthesize individual files into group files and remove missing variants.
# 
# Input:
# - VCF files, genetic maps, and population-specific configuration files.
# 
# Output:
# - Local ancestry results, integrated ancestral haplotypes, and assembled VCF files.

#Software
python=/your/path			#Python3
plink1.9=/your/path			#Plink v1.9
plink2.0=/your/path			#Plink v2.0
perl=/your/path				#Perl v5
ChromoPainterv2=/your/path	#ChromoPainter v2
bcftools=/your/path			#BCFtools v1.14


echo "Step 1: Local ancestry inference for the target populations"
# ChromoPainterV2 uses Hidden Markov Model analysis to intuitively form each recipient chromosome as a mosaic of the donor chromosomes, capturing which donors are essential to explain the recipient DNA

#(1). Make ChromoPainterV2 input file of genetic maps. It contains the recombination distance between each pair of contiguous SNVs
#Process genetic maps to the input format of Plink (v1.9)
zcat   chr21.b38.gmap.gz   |   sed   '1d'   |   awk   '{print $2"\t"$1"\t"$3}'   >   101.shapeit4.chr21.b38.map
#python   01.recombination.py   <Input>   <Output>
#<Input> is a tab-separated file with no header. Column 1: Start position of genomic markers; Column 2: Genomic marker identifier; Column 3: Genetic map distance (in cM).
#<Output> is a space-separated file containing: Start position; Recombination rate per base pair; Recombination rate in cM per megabase; Genetic map distance
python   01.recombination.py   101.shapeit4.chr21.b38.map   102.shapeit4.chr21.b38.recomb.txt
echo   "position . COMBINED_rate(cM/Mb) Genetic_Map(cM)"   |   cat   -   102.shapeit4.chr21.b38.recomb.txt   |   sed '2d'   |   awk  '{print $1,$3,$4}'   >   103.genetic_map_chr21_combined_b38.txt

#Add genetic distance for each SNV in test dataset
plink1.9   --vcf Test.dataset.chr21.vcf.gz   --double-id   --cm-map 103.genetic_map_chr21_combined_b38.txt   21   --make-bed   --out 111.chr21.samples.filter.phased

#Make ChromoPainterV2 input (recombination rate infile)
awk   '{print $1"\t"$4"\t"$3}'   111.chr21.samples.filter.phased.bim   >   121.chr21.samples.filter.phased.input.txt
python   01.recombination.py   121.chr21.samples.filter.phased.input.txt   122.chr21.samples.filter.phased.recomb.txt
awk   '{print $1" "$2}'   122.chr21.samples.filter.phased.recomb.txt  >   123.chr21.samples.filter.phased.recomrates

#(2). Make ChromoPainterV2 input file of haplotypes. It contains the SNV data for a set of chromosomes. For 'impute2chromopainter2.pl', The '-r' parameter of the script needs to be followed by the number of samples in the input file. Our test data has 25 samples.
plink2.0   --vcf Test.dataset.chr21.vcf.gz   --double-id   --export haps ref-first   --chr 21   --out 131.chr21.samples.filter.phased
perl   02.impute2chromopainter2.pl   -r 25   131.chr21.samples.filter.phased.haps   132.chr21.samples.filter.phased

#(3). Make ChromoPainterV2 input file of individual labels. The default use of all samples in the test dataset.
awk   '{print $1,$3,"1"}'   sample.list.txt   >   141.idfile.txt

#(4). Make ChromoPainterV2 input file of population lists. It contains recipient and donor populations. As an example, we choose to reconstruct the ancestral genomes of the Miao and She ethnic groups, which we refer to as target populations. We have prepared two files of population lists. In the following two sets of information, these two files are for recipient and donor populations.
echo 'Miaozu D' > 151.poplist.txt
echo 'She R' >> 151.poplist.txt
echo 'Han_NChina D' >> 151.poplist.txt
echo 'Han D' >> 151.poplist.txt
echo 'Dai D' >> 151.poplist.txt

echo 'Miaozu R' > 152.poplist.txt
echo 'She D' >> 152.poplist.txt
echo 'Han_NChina D' >> 152.poplist.txt
echo 'Han D' >> 152.poplist.txt
echo 'Dai D' >> 152.poplist.txt

#(5). Run ChromoPainterV2 to perform local ancestry inference. For the two population files, run ChromoPainterV2 separately. By default, haplotypes for each recipient are inferred 10 times. To save storage, compress the output files '*.samples.out'
ChromoPainterv2   -g 132.chr21.samples.filter.phased.haps \
		-r 123.chr21.samples.filter.phased.recomrates \
		-t 141.idfile.txt \
		-f 151.poplist.txt 0 0 \
		-o 161.cp.chr21 -i 10 -in -iM \
		&& bgzip   -f   -@ 5   161.cp.chr21.samples.out

ChromoPainterv2   -g 132.chr21.samples.filter.phased.haps \
		-r 123.chr21.samples.filter.phased.recomrates \
		-t 141.idfile.txt \
		-f 152.poplist.txt 0 0 \
		-o 162.cp.chr21 -i 10 -in -iM \
		&& bgzip   -f   -@ 5   162.cp.chr21.samples.out




echo "Step 2: Determine which ancestries the recipient haplotypes copies from at each allele"
# For the 10 inferences, if a donor population is identified at least 5 times, the ancestry of the recipient haplotype at that variant is assigned to the donor population

# Create ID file for the analysis
awk   'BEGIN{R=1}{print $1"\t"$2"\t"R; R+=1; print $1"\t"$2"\t"R; R+=1}'   141.idfile.txt   >   201.id.txt

# "03.identify.ancestry.py" determines the ancestry of haplotypes using ChromoPainter output.
# python   03.identify.ancestry.py   <ChromoPainter output file>  <Population information file>   <Output filename>   <Cutoff>
num=5
python   03.identify.ancestry.py \
	161.cp.chr21.samples.out.gz \
	201.id.txt \
	211.chr21.part1.${num}.txt.gz \
	${num}

python   03.identify.ancestry.py \
	162.cp.chr21.samples.out.gz \
	201.id.txt \
	211.chr21.part2.${num}.txt.gz \
	${num}



echo "Step 3: Data transpose and adding physical location information"
# Data transpose and adding physical location information

num=5
# "04.transpose.py" converts the ancestry information file from a position-centric format to a sample-centric format.
# python   04.transpose.py   <Ancestry information with positions as columns>   <position list file>   <Ancestry information with samples as columns>
python   04.transpose.py \
	211.chr21.part1.${num}.txt.gz \
	131.chr21.samples.filter.phased.haps \
	221.chr21.ind.pos.part1.${num}.txt.gz

python   04.transpose.py \
	211.chr21.part2.${num}.txt.gz \
	131.chr21.samples.filter.phased.haps \
	221.chr21.ind.pos.part2.${num}.txt.gz




echo "Step 4: Extract the specified ancestral alleles from the selected samples"
# Extract the specified ancestral alleles from the selected samples

# "05.extract.identified.ancestral.alleles.py" extracts alleles belonging to the specified ancestry from the selected samples in phased VCF files.
# python   05.extract.identified.ancestral.alleles.py   <Ancestry information with samples as columns>   <selected samples>   <Phased VCF>   <Output ancestral gene pool>   <Specified ancestry to be extracted>
# Separately extract alleles identified as belonging to the Miao ancestry in She samples and alleles identified as belonging to the She ancestry in Miao samples.
python   05.extract.identified.ancestral.alleles.py   221.chr21.ind.pos.part1.5.txt.gz \
	sample.list.She.txt \
	Test.dataset.chr21.vcf.gz \
	231.chr21.5.She.txt.gz \
	Miaozu

python   05.extract.identified.ancestral.alleles.py   221.chr21.ind.pos.part2.5.txt.gz \
	sample.list.Miaozu.txt \
	Test.dataset.chr21.vcf.gz \
	231.chr21.5.Miaozu.txt.gz \
	She


echo "Step 5: Assemble ancestral haplotypes using extracted ancestral gene pools"
# Assemble ancestral haplotypes using extracted ancestral gene pools
cutoff_freq=0.2
cutoff_length=6000
num=5

# "06.inte.genome.cutoff.py" integrates ancestral gene pool across multiple files to generate an ancestral haplotype.
# python   06.inte.genome.cutoff.py   <Output ancestral haplotype>   <Missing rate cutoff>   <Length cutoff>   <Ancestral gene pool 1>   …   < Ancestral gene pool n>
# Length cutoff: The maximum length that can be extended at one time when obtaining ancestral segments on a haplotype (e.g., 6000.0 bp).
# Can input one or more ancestral gene pool files by appending them to the script parameters in sequence.

for i in {1..10}    # Generate 10 ancestral haplotypes
do
	python   06.inte.genome.cutoff.py  \
	301.Miao.She.chr21.${i}.${cutoff_freq}.${cutoff_length}.${num}.txt.gz \
	${cutoff_freq} \
	${cutoff_length} \
	231.chr21.5.Miaozu.txt.gz \
	231.chr21.5.She.txt.gz
done



echo "Step 6: Generate ancestral individuals using assembled ancestral haplotypes"
# Generate ancestral individuals from the assembled haplotypes
num=5
cutoff_freq=0.2
cutoff_length=6000

# "07.trans.to.vcf.py" converts two ancestral haplotypes into a VCF-format ancestral individual.
# python   07.trans.to.vcf.py   <Input gzipped VCF file>   <Ancestral haplotype 1>   < Ancestral haplotype 2>   <Sample ID>   <Output ancestral individual>

for i in `seq 1 2 10`
do
	j=$[ $i + 1 ]
	#Combine two consecutively numbered haplotypes to generate a diploid
	python   07.trans.to.vcf.py \
	Test.dataset.chr21.vcf.gz \
	301.Miao.She.chr21.${i}.${cutoff_freq}.${cutoff_length}.${num}.txt.gz \
	301.Miao.She.chr21.${j}.${cutoff_freq}.${cutoff_length}.${num}.txt.gz \
	MiaoShe.${i}.${j} \
	401.MiaoShe.chr21.${i}.${j}.${cutoff_freq}.${cutoff_length}.${num}.vcf
done



echo "Step 7: Synthesize individual files into group files and remove missing variants"
# Synthesize individual files into group files and remove missing variants
cutoff_freq=0.2
cutoff_length=6000
num=5

# Loop through chromosomes to merge individual VCF files
for   chr   in   21
do
	rm   501.merge.chr${chr}.list.txt
	
	# Generate a list of VCF files for merging
	for  i   in   `seq 1 2 10`
	do
		j=$[ i + 1 ]
	
		echo   "401.MiaoShe.chr${chr}.${i}.${j}.${cutoff_freq}.${cutoff_length}.${num}.vcf.gz"   >>   501.merge.chr${chr}.list.txt
	done

	# Merge the VCF files
	bcftools   merge   -l   501.merge.chr${chr}.list.txt   -O z   -o 502.assemble.chr${chr}.${cutoff_freq}.${cutoff_length}.${num}.tmp.vcf.gz   &&    tabix   -f   -p vcf   502.assemble.chr${chr}.${cutoff_freq}.${cutoff_length}.${num}.tmp.vcf.gz

	# Rename chromosomes in the merged file.
	echo   -e   "${chr}\tchr${chr}"   >   503.rename_chr${chr}.txt
	bcftools   annotate   --rename-chrs   503.rename_chr${chr}.txt   -O z   -o 504.assemble.chr${chr}.${cutoff_freq}.${cutoff_length}.${num}.vcf.gz   502.assemble.chr${chr}.${cutoff_freq}.${cutoff_length}.${num}.tmp.vcf.gz   &&   tabix   -f   -p vcf   504.assemble.chr${chr}.${cutoff_freq}.${cutoff_length}.${num}.vcf.gz

	# Remove missing variants
	bcftools   view   -i 'F_MISSING<0.001'   --threads 10   -O z   -o 505.assemble.chr${chr}.${cutoff_freq}.${cutoff_length}.${num}.rm.missing.vcf.gz  504.assemble.chr${chr}.${cutoff_freq}.${cutoff_length}.${num}.vcf.gz   &&   tabix   -f   -p vcf   505.assemble.chr${chr}.${cutoff_freq}.${cutoff_length}.${num}.rm.missing.vcf.gz
done

