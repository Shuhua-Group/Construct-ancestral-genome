#!/usr/bin/env Rscript
#### NEED PACKAGES 'data.table', 'Rcpp', 'stringr', optparse' FOR THIS PROGRAM TO WORK ####
#### IF RUNNING IN UNCERTAINTY MODE, THERE NEEDS TO BE 2 DIFFERENT VCF FILES - ONE OF THEM CONTAINING 
#### If you want to print a recomrates file, you must have the "CM" field in the INFO tag of the genotypes file - SHAPEIT4 does this automaticall the tag must also be exactly named  CM=x and the field seperated by ";" ##

library(data.table)
library(stringr)

option_list = list(
	optparse::make_option(c("-g", "--genotypes"), 
		type="character", 
		default=NULL, 
		help="genotype file name. no default - required", 
		metavar="character"), 
	optparse::make_option(c("-l", "--genotypelikelihoods"), 
		type="character", 
		help="genotype likelihoods file name", 
		metavar="character"), 
	optparse::make_option(c("-u", "--uncertaintyMode"), 
		type="logical", 
		default="F", 
		help="Run in uncertainty mode", 
		metavar="logical"), 
	optparse::make_option(c("-o", "--output"), 
		type="character", 
		default=NULL, 
		help="Output file stem. no default - required", 
		metavar="character"),
	optparse::make_option(c("-m", "--genmap"),
		type="character",
		default=NULL,
		help="plink formatted genetic map",
		metavar="character")
    );


source_cpp_f = function() {
  Rcpp::sourceCpp("/home/zhangxiaoxi/software/vcf_to_chromopainter-main/src/vcf_to_chromopainter_functions.cpp")
  cat("Successfully souced and compiled source code!\n")
}

make_recombination_map = function(geneticmap, positions) {
	genmap = fread(geneticmap)
	positions = sort(positions)
	data.table(start.pos = positions, recom.rate.perbp = approx(x = genmap$V4, y = genmap$V3, xout = positions)$y)
}

readVCF = function(likelihoods_file) {
	fread(cmd = paste("zgrep -v '^##'", likelihoods_file))
}

getDSField = function(formatField){
	str_which(str_split_1(formatField, ":"), "DS")
}

getGTField = function(formatField){
	str_which(str_split_1(formatField, ":"), "GT")
}

getGPField = function(formatField){
	str_which(str_split_1(formatField, ":"), "GP")
}

processGenotypes = function(filename) {
	genotypes = readVCF(filename)
	POS = genotypes$POS
	if (colnames(genotypes)[1] != "#CHROM") {
		stop("First field isn't chromosome - skipped the wrong number of lines? Exiting....")
	}
	subsetStart = which(colnames(genotypes) == "FORMAT")
	dropcols = colnames(genotypes)[1:subsetStart]
	list(POS, genotypes[, (dropcols) := NULL])
}

processGenotypeLikelihoods = function(filename) {
	likelihoods = readVCF(filename)
	if (colnames(likelihoods)[1] != "#CHROM") {
		stop("First field isn't chromosome - skipped the wrong number of lines? Exiting....")
	}
	field = str_which(likelihoods[1,], "GP")
	format = as.character(likelihoods[1,..field])
	dsField = getDSField(format)
	gtField = getGTField(format)
	gpField = getGPField(format)
	GLpositions = likelihoods$POS
	subsetStart = which(colnames(likelihoods) == "FORMAT")
	dropcols = colnames(likelihoods)[1:subsetStart]
	likelihoods = likelihoods[, (dropcols) := NULL]
	field1 = likelihoods[1,1]
	gpField1 = str_split_1(as.character(field1), ":")[gpField]
	gpdelim = gpField1 %>% str_remove_all("\\.") %>% str_extract("[:punct:]")
	list(likelihoods, gpdelim, gpField-1)
}


###################### Do stuff ##############################################

opt_parser = optparse::OptionParser(option_list=option_list);
opt = optparse::parse_args(opt_parser);

if (is.null(opt$output)) {
	optparse::print_help(opt_parser)
	stop("Must supply output file name\n", call.=FALSE)
	}

if (is.null(opt$genotypes)) {
	optparse::print_help(opt_parser)
		stop("At least one argument must be supplied (input file)\n", call.=FALSE)
	}

if (opt$uncertaintyMode == TRUE) {
	cat("Running in Uncertainty mode!\n")
}

if (!opt$uncertaintyMode && is.null(opt$genotypelikelihoods)) {
  mode = "first"
	cat("Running in classic chromopainter mode!\n")
} else if (opt$uncertaintyMode && is.null(opt$genotypelikelihoods)) {
  mode = "second"
	cat("Running in uncertainty mode but not accounting for genotype likelihoods\n")
} else if (opt$uncertaintyMode && !is.null(opt$genotypelikelihoods)) {
  mode = "third"
	cat("Running in uncertainty mode and also accounting for genotype likelihoods\n")
} else {
	stop("Combination of parameters is incorrect")
}


CPoutname = paste0(opt$output, ".chromopainter.inp")
recomratesoutname = paste0(opt$output, ".recomrates.txt")
idfileout = paste0(opt$output, ".idfile.txt")

source_cpp_f()

GTs = processGenotypes(opt$genotypes)
POS = GTs[[1]]
GTs = GTs[[2]]

cat("Making recombination rates file\n")
recomap = make_recombination_map(opt$genmap, POS)
fwrite(recomap, recomratesoutname, sep=" ")
cat("Finished making recombination rates file\n")

cat("Making idfile\n")
idfile = data.table(V1 = colnames(GTs), V2 = colnames(GTs), V3 = 1)
fwrite(idfile, idfileout, sep=" ")
cat("Finished making idfile\n")

### typicical uncertainty mode 
if (mode == "third") {
	GLs = processGenotypeLikelihoods(opt$genotypelikelihoods)
	gpdelim = GLs[[2]]
	GPfield = GLs[[3]]
	GLs = GLs[[1]]
	CPout = ReturnChromopainterUncerainty(as.matrix(GTs), as.matrix(GLs), GPfield, gpdelim)
	POS[1] = paste0("P ", POS[1])
	colnames(CPout) = POS
	fwrite(CPout, CPoutname, sep=" ", quote=F)
}

## this is if we don't have genotypes but we need to print them out in uncertainty format ##
## then we write out in uncertainty format but with no actual uncertainty in the values ##
if (mode == "second") {
	CPout = ReturnChromopainter(as.matrix(GTs))
	CPout = as.data.table(CPout)	
	POS[1] = paste0("P ", POS[1])
	colnames(CPout) = POS
	fwrite(CPout, CPoutname, sep=" ", quote=F)	
}

## just typical chromopainter with genotypes only and no uncertainty ##
if (mode == "first") {
	CPout = ReturnChromopainter(as.matrix(GTs))
	CPout = as.data.table(CPout)
	CPout = CPout[, combined := do.call(paste, c(.SD, sep = ""))][, "combined", with=F]
	POS[1] = paste0("P ", POS[1])
	CPout = rbindlist(list(data.table(combined = paste(POS, collapse=" ")), CPout))
	fwrite(CPout, cpoutname, sep = "\n", col.names = FALSE, quote = FALSE)
}

## write out the annoying chromopainter header lines ##
CPoutfile = readLines(CPoutname)
headlines = paste0(length(POS), "\n", ncol(GTs) * 2)
CPoutfile = c(headlines, CPoutfile)
writeLines(CPoutfile, CPoutname)

cat("Completed!\n")
