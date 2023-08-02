##### SKIP THIS STEP IF YOU HAVE ALREADY DOWNLOADED THESE PACKAGES BEFORE THE COURSE

#install.packages("ggplot2")
#install.packages("dada2")
#install.packages("ShortRead")
#install.packages("Biostrings")
#install.packages("Rcpp")

# Load packages
library(ggplot2)
library(dada2)
library(ShortRead)
library(Biostrings)


# Define path to files

path <- "~/Documents/Leticia_data_2"  ## CHANGE ME to the directory containing the fastq files.
list.files(path)

# Generate lists of forward and reverse files 

fnFs <- sort(list.files(path, pattern = "R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "R2_001.fastq.gz", full.names = TRUE))


## Plot quality profiles of the Forward reads

plotQualityProfile(fnFs[1:6])

## Plot quality profiles of the Reverse reads

plotQualityProfile(fnRs[1:6])


### Input primer sequences (ITS3 and ITS4) ###

FWD <- "GCATCGATGAAGAACGCAGC"  ## CHANGE ME to your forward primer sequence
REV <- "TCCTCCGCTTATTGATATGC"  ## CHANGE ME...

## Check these are the correct primers (and that they are found in the correct orientation in F and R)

allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients


# Filter to remove Ns

fnFs.filtN <- file.path(path, "filtN", basename(fnFs)) # Put N-filtered files in filtN/ subdirectory
fnRs.filtN <- file.path(path, "filtN", basename(fnRs))
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE)



### Count primer occurrences ###

primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]]))


# Cutadapt installed with conda (conda activate mambaforge/envs/cutadaptenv )
cutadapt <- "/Users/bethanmanley/mambaforge/bin/cutadapt"
system2(cutadapt, args = "--version") # Run shell commands from R


# Run cutadapt

path.cut <- file.path(path, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))

FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD, "-a", REV.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV, "-A", FWD.RC) 
# Run Cutadapt
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                             fnFs.filtN[i], fnRs.filtN[i])) # input files
}


# Sanity check - count presence of primers in first sample post cutadapt

rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]]))


##########################################################################

## Analyse primer-free sequences through cutadapt

# Forward and reverse fastq filenames have the format:
cutFs <- sort(list.files(path.cut, pattern = "_R1_001.fastq.gz", full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern = "_R2_001.fastq.gz", full.names = TRUE))

# Extract sample names, assuming filenames have format:
get.sample.name <- function(fname) strsplit(basename(fname), "_S")[[1]][1]
sample.names <- unname(sapply(cutFs, get.sample.name))
head(sample.names)



## Assigning output names

filtFs <- file.path(path.cut, "filtered", basename(cutFs))
filtRs <- file.path(path.cut, "filtered", basename(cutRs))


## Filtering

out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, maxN = 0, maxEE = c(2, 4), 
                     truncQ = 2, minLen = 50, rm.phix = TRUE, compress = TRUE, multithread = TRUE)  # on windows, set multithread = FALSE
head(out)


# Inspect read quality profiles of forward and reverse reads

plotQualityProfile(filtFs[1:6])

plotQualityProfile(filtRs[1:6])



#### Error rate estimation using the data ####

errF <- learnErrors(filtFs, multithread=TRUE)

errR <- learnErrors(filtRs, multithread=TRUE)


# Visualise estimated error rates

plotErrors(errF, nominalQ=TRUE)


################ Sample inference ################

# Applies the core sample inference algorithm to the filters/trimmed seqs
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)

dadaRs <- dada(filtRs, err=errR, multithread=TRUE)

# Check the returned object (dada-class)
dadaFs[[8]]


## By default, DADA2 processes and infers sequence varieties from each sample independently, but you can choose to pool samples to detect rarer sequence variants.

######### Merge paired reads ########

# Align denoised forward and reverse reads, construct merged contig sequences. Default paraeter is exact matches of 12bp overlap between F and R

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

# Inspect the merger data.frame from the first sample
head(mergers[[8]])


#### Construct ASV table ######

seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))


## Remove chimeras ##

# Identified if they can be exactly reconstructed by combining a left-segment and a right-segment from two more abundant "parent" sequences

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

sum(seqtab.nochim)/sum(seqtab)


## Track number of reads that made it through each step

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
track


#### Assign taxonomy ####

## Requires a set of input sequences to be classified and a training set of reference sequences with known taxonomy. For fungi, use the General Fasta release files from the UNITE ITS database -https://unite.ut.ee/repository.php#general. 

taxa <- assignTaxonomy(seqtab.nochim, "~/Documents/Data/UNITE_general_release/sh_general_release_dynamic_18.07.2023.fasta", multithread=TRUE)


## Inspect taxonomic assignments

taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)


### There is an alternative to the assignTaxonomy function (described here: https://benjjneb.github.io/dada2/tutorial.html) - apparently assigns taxonomy faster and possibily better



####### Handoff to phyloseq ########

library(phyloseq); packageVersion("phyloseq")

library(Biostrings); packageVersion("Biostrings")

library(ggplot2); packageVersion("ggplot2")

theme_set(theme_bw())


### Read in sample metadata from file

samdf <- read.csv("~/Documents/Leticia_data_2/leticia_vegetation.csv")

samples.out <- rownames(seqtab.nochim)
rownames(samdf) <- samples.out

## Construct phyloseq object

physeq_dada2 <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa))


### Store DNA sequence in refseq slot, rename taxa to a short string

dna <- Biostrings::DNAStringSet(taxa_names(physeq_dada2))
names(dna) <- taxa_names(physeq_dada2)
physeq_dada2 <- merge_phyloseq(physeq_dada2, dna)
taxa_names(physeq_dada2) <- paste0("ASV", seq(ntaxa(physeq_dada2)))
physeq_dada2
