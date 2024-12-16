# Comparative Genomics ST1 *Staphylococcus aureus*
![bash](https://img.shields.io/badge/language-bash-green)
![Python](https://img.shields.io/badge/language-Python-blue)
![R](https://img.shields.io/badge/language-R-red)

All bioinformatic analysis was conducted on the New Zealand eScience Infrastructure [NeSI](https://github.com/nesi). FastQC and Kraken2 was used for QC of the samples. A hybrid genome assembly was used to create the first ST1 bovine *S. aureus* reference genome 23EV612. Nullarbor was used for genome analysis, SpaTyper and AGRvate was used to determine *spa* and *agr* types. ST1 Phylogenetic tree generation using Gubbins and IQtree2. Python was used to calculate Inter and Intra-host variation between bovine and human *S. aureus* genomes and also calculate statistically significant genes based on presence and absence of detected antimicrobial resistance genes, Virulence and mobile genetic elements (MGE). Characterisation and determination of ΦSabovST1 is documented below. 

## FastQC 
[FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) was used to check for QC of the samples for adaptor content and sequence quality
```#!/bin/bash -e
#SBATCH --cpus-per-task=8 --mem 50Gb --time 1:00:00 -J FASTQC_EV

module load FastQC/0.11.9

FASTQC = /pathtorawreads
OUTPUT_DIR = /processedreadsdirectory/FASTQC/

mkdir -p $OUTPUT_DIR

for file in $FASTQC/*.fq.gz;
do
if [-f $FILE ];
then echo "Running FastQC on $FILE"
fastqc $FILE -o $OUTPUTDIR
else
echo "No FASTQ files found in $FASTQ"
fi
done
echo "FastQC analysis completed for all samples"
```

## Kraken2 
The default Kraken2 installed with [Nullarbor](https://github.com/tseemann/nullarbor) was used for species identification. If an isolate had a *k-mer* match ≤70% to *S. aureus* and/or contained more than 5% *k-mer* matches to another bacterial species these samples will be excluded. 

```#!/bin/bash
CONF=0.5
DATA=`find ${1} -name "*R1*.fq.gz"`
mkdir -p Kraken_Reports/
OUTPUT= Kraken_Reports/

for i in ${DAT[@]}
do
        FQ1=${i}
        FQ2=`echo ${i} | sed s/R1/R2/`
        echo $FQ1
        echo $FQ2
        base=$(basename ${i}  _R1.fq.gz)
        echo "basename is ${base}"
        echo "#!/bin/bash" > tmp.sh
        echo "set -x; module load nullarbor/2.0.20191013; kraken2 --confidence ${CONF} --report ${OUTPUT}/${base}.report --threads 24 --paired ${FQ1} ${FQ2} > ${OUTPUT}/${base}.log" >> tmp.sh
        sbatch -J KrakenIllumina_EV --mem 50G --time 0-1 -c 24 tmp.sh
        sleep 0.5
done
```

## Assembly of 23EV612 
The Hybrid Genome Assembly of *S. aureus* isolate 23EV612 is detailed below. 
#### Guppy
[Guppy](https://community.nanoporetech.com/docs/prepare/library_prep_protocols/Guppy-protocol/v/gpb_2003_v1_revax_14dec2018/guppy-software-overview) was chosen as the basecaller specifying the super accurate model. 
#### Filtlong
[Filtlong](https://github.com/rrwick/Filtlong) was used to remove reads shorter than 1kbp and exclude the worst 5% of the reads. 
#### Porechop
[Porechop](https://github.com/rrwick/Porechop) was used to remove any adaptors 
#### NanoStat
[NanoStat](https://github.com/wdecoster/nanostat) was used to perform quality checks to determine the mean read quality of 23EV612
NanoStat Output 
|General Summary  |Quality Statistics |
|------------------|--------------------|
|Mean Read Length  | 6053.4             |
|Mean Read Quality | 17.8               |
|Median Read Length | 5741.0            | 
|Median Read Quality | 17.5             | 
|Number of Reads  | 1,262,085.0        |
|Number of bases  | 7,689,967,881      |
|Read length N50 | 7,549.0             | 
|>Q5             | 100% 7640.0Mb      |
|>Q7             | 100% 7640.0Mb      |
|>Q10            | 100% 7640.0Mb      |
|>Q12            | 90.7% 6866.0Mb      |
|>Q15            | 58.9% 4450.1Mb      |

#### Chopper
[Chopper](https://github.com/wdecoster/chopper) was used to ensure only reads above the mean read quality are kept. 
#### Flye 
[Flye](https://github.com/mikolmogorov/Flye) was used as the *de novo* assembler to assemble the isolate 23EV612. 
#### Medaka 
[Medaka](https://github.com/nanoporetech/medaka) was used as the polishing tool for 23EV612. 
#### BWA and SAMtools 
[BWA](https://github.com/lh3/bwa) and [SAMtools](https://samtools.sourceforge.net/) were used to iteratively map the Illumina reads of 23EV612 to the nanopore assembly. 
#### Pilon
[Pilon](https://github.com/broadinstitute/pilon) is run to polish the final Hybrid Assembly 23EV612. Three rounds of Pilon are run. Therefore Pilon is run three times and the bwamem_pilon is run twice. The final fasta file called pilon_round3.fasta is the polished hybrid genome and is renamed accordingly. 
#### CheckM and Quast 
[CheckM](https://github.com/Ecogenomics/CheckM) and [Quast](https://github.com/ablab/quast) was run on the 23EV612 to check the assembly is complete 

## Nullarbor
[Nullarbor](https://github.com/tseemann/nullarbor) is the genomic analysis pipeline chosen for bovine, human and ST1 *S. aureus* analysis. In the nullarbor.pl file --mincov was specified as 50 and --minid as 90. The script was run on NESI. 
```#!/bin/bash -e
#SBATCH --cpus-per-task=20 --mem 160Gb --time 166:00:00 -J nullarbor_EV
#Alter Reference Genome for human, bovine or ST1 analysis - comment out the one not needed and alter Reference in the script.
#Alter text file for purpose - bovine, human or ST1 analysis
##Ref_Bovine = SaureusRF122.fa
#Ref_Human = SaureusMSSA476.fa
Ref_ST1 = 23EV612.fa

mkdir -p $TMPDIR
module purge
module load nullarbor/2.0.20191013

nullarbor.pl --name StaphAureusIsolates --ref $Ref_ST1 \
--input ST1_Isolates.txt \
--outdir ST1_Analysis --force \
--cpus $SLURM_CPUS_PER_TASK --run --mlst saureus --trim --taxoner-opt confidence 0.5 --treebuilder iqtree --treebuilder-opt "-wbtl -st DNA -b 1000 --alrt 1000" --annotator-opt "--fullanno --compliant"
```

## Defining *Spa* Types 
*Spa* Types were determined using [SpaTyper](https://github.com/HCGB-IGTP/spaTyper). The list of genomes is made from the contigs.fa created by Nullarbor. 
```ls *.fa > list_of_genomes.txt ##Creating a list of genomes from all fastq files 
sed -i '1 i\spaTyper  -f ' list_of_genomes.txt
echo "--output spa.xls" &>> list_of_genomes.txt
tr '\n' ' ' < list_of_genomes.txt > spa_command.txt
chmod +x spa_command.txt
./ spa_command.txt
```
## Defining *agr* types
[AGRvate](https://github.com/VishnuRaghuram94/AgrVATE) was used to identify each genomes *agr* type. A conda environment was created then the contig fasta files was used as input for AGRvate. 
## Gubbins
[Gubbins](https://github.com/nickjcroucher/gubbins) was used to remove recombinant regions from the ST1 alignment file before constructing a time-scaled tree. 
The alignment file from Nullarbor needed to be 'cleaned' using the snippy-clean_full_aln command installed with [Snippy](https://github.com/tseemann/snippy). Gubbins script was run as detailed below. 
```snippy-clean_full_aln core.full.aln > clean_full.aln```
```#!/bin/bash -e
#SBATCH --cpus-per-task=8 --mem 40Gb --time 166:00:00 -J Gubbins_EV -e Gubbins_%J.err -o Gubbins_%J.out

mkdir -p $TMPDIR
module purge
module load Gubbins/3.2.2-gimkl-2022a-Python-3.10.5
module load RAxML-NG/1.1.0-gimkl-2022a

run_gubbins.py -v --prefix ST1_NZ_Gubbins clean.full.aln  --first-tree-builder rapidnj --first-model JC -# 1000 
```
## IQ-TREE2
[IQ-TREE2](https://github.com/iqtree/iqtree2) was chosen to generate a time-scaled analysis ST1 NZ tree using the recombination sites and tree as output from Gubbins as input to IQ-TREE2 with the addition of a date file which is tab-delimited - isolate name in column 1 and date in column 2 enusring dates in YYYY-MM-DD format. 
```#!/bin/bash -e
#SBATCH --cpus-per-task=8 --mem 50Gb --time 166:00:00 -J IQTREE_EV
module load nullarbor/2.0.201913

iqtree2 --date ST1_NZ_dates.txt -s CC1_TimeTree_Gubbins.filtered_polymorphic_sites.phylip --tree CC1_TimeTree_Gubbins.filtered_polymorphic_sites.phylip.treefile
```
## Calculating Inter and Intra-host SNP variation using the SNP distance matrix 
The SNP distance matrix as output from Nullarbor was used along with a CSV file which contained isolates in column 1 and column 2 was the respective host. Python code was used to compute the SNP distances between human and bovine hosts. Script called InterIntraHost.py

## Calculating Statistical Significance Antimicrobial Resistance and Virulence genes and also Mobile Genetic Elements 
Python script was used to calculate if genes were statistically significant (StatisticallySignficanceGenes.py) using Chi-Square/Fisher's Exact based on presence/absence data of each detected gene in human and bovine hosts. 

## Virulence Gene Analysis 
Virulence Genes were split into their broad function: Adherence, Enterotoxin, Exoenzyme, Exotoxin, Haemolysin, Immune Modulation, Intracellular Adhesion, Type VII Secretion System and Others. The csv file was a binary matrix where Isolate in column 1, host in column 2 and the following columns with the genes detected where 1 is present and 0 being absent. Each individual function csv file was the input into the below python script to put *S. aureus* isolates of each host into groups to identify the combination of genes that are present in majority of each bovine and human *S. aureus* host. Each functional group output was visualised as a presence/absence heatmap in R. 

```python
file_path = "exoenzyme.csv" #Altered for each functional group
results = {}

with open(file_path,'r') as data:
    for line in data.readlines():
        if 'Host' in line:
            continue
        else:
            tmpdat = line.split(',')
            key = line[line.index(',') + 1: len(line)].rstrip()
            if key in results:
                results[key].append(tmpdat[0])
            else:
                results[key] = [tmpdat[0]]

print()
for items in results:
    print("{},{},{}".format(items,len(results[items]),','.join(results[items])))
```
## Detection of φSabovST1
The phage sequence φSabovST1 was determined as a similar match to φSaov3 in our MGE database built using abricate. In the abricate output, genome coordinates for the match are displayed, these genome coordinates were used to extract the sequence of interest. However, 10,000bp either side of the detected region was chosen as a flanking region. The phage region was uploaded to [Phastest](https://phastest.ca/) and [PhageScope](https://phagescope.deepomics.org/) to determine the structure of the detected phage. Core-SNP phylogenetic tree between the reference φPV83, φSaov3, φSabovST1 and [Cluster B7 Phages](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-019-5647-8) was determined using [snippy](https://github.com/tseemann/snippy), [snp-dists](https://github.com/tseemann/snp-dists) and [IQ-TREE2](https://github.com/iqtree/iqtree2). [FastANI](https://github.com/ParBLiSS/FastANI) was used to determine the mean average nucleotide identity for all the phages used in the phylogeny. The ANI output was visualised in R.  
```module load SAMtools/1.16.1-GCC-11.3.0
samtools faidx 23EV612.fasta Contig1:2730000-2780000 > 23EV612_phagesequence.fasta
```
```#!/bin/bash -e
#SBATCH --cpus-per-task=20 --mem 50Gb --time 1:00:00 -J SNIPPY_EV
mkdir -p $TMPDIR
module purge
module load nullarbor/2.0.20191013

REF=phiPV83.fa
INPUT=Phage.tab #Tab delimited file with phage in column 1 and path to sequence in column 2. 

snippy-multi $INPUT --ref $REF --cpus 16 > snippy_phages.sh
echo snippy_phages.sh has been made
sh ./snippy_phages.sh
```
```snippy-clean_full_aln core.full.aln > clean_full.aln```
```snp-dists clean_full.aln > distances.tab```
```#!/bin/bash -e
#SBATCH --cpus-per-task=8 --mem 50Gb --time 1:00:00 -J IQTREE_EV
module load nullarbor/2.0.201913
iqtree2 -s clean_full.aln 
```
```#!/bin/bash
#SBATCH -J FastANI
#SBATCH --time 00:30:00
#SBATCH --mem 10GB
#SBATCH --ntasks=1
#SBATCH --array=0-13
#SBATCH --cpus-per-task=15

# Working directory
cd FASTANI/

# Load module
module purge
module load FastANI/1.33-GCC-11.3.0

#Creating the array
MY_FILE_LIST=Phage.txt
LIST_OF_ALL_FILES=(`cat ${MY_FILE_LIST}`)
ARRAY_ITEM=${LIST_OF_ALL_FILES[$SLURM_ARRAY_TASK_ID]}

# Run FastANI
echo ${ARRAY_ITEM}
fastANI -r ${ARRAY_ITEM} --ql $MY_FILE_LIST -o ${SLURM_ARRAY_TASK_ID}.txt
```








