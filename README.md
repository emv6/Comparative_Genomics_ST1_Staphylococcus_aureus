# Comparative Genomics ST1 *Staphylococcus aureus*

All bioinformatic analysis was conducted on the New Zealand eScience Infrastructure [NeSI](https://github.com/nesi). FastQC and Kraken2 was used for QC of the samples. A hybrid genome assembly was used to create the first ST1 bovine *S. aureus* reference genome 23EV612. Nullarbor was used for genome analysis, SpaTyper and AGRvate was used to determine *spa* and *agr* types. ST1 Phylogenetic tree generation using Gubbins and IQtree2. Python was used to calculate Inter and Intra-host variation between bovine and human *S. aureus* genomes and also calculate statistically significant genes based on presence and absence of detected antimicrobial resistance genes, Virulence and mobile genetic elements (MGE). Characterisation and determination of ΦSabovST1 is documented below. 

## FastQC 
[FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) was used to check for QC of the samples for adaptor content and sequence quality
```bash
#!/bin/bash -e
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

```bash
#!/bin/bash

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
```bash 
#!/bin/bash -e
#SBATCH -J guppy_gpu --gpus-per-node A100:1 --mem 6G --cpus-per-task 4 --time 10:00:00 --output slurmout.%j.out

module purge
module load ont-guppy-gpu/6.4.6

INPUT=Path/to/raw/reads/20201219_1143_MC-110462_0_FAN52032_782f79b0/
FOLDER=Guppy_Super_Accurate_model_23EV612/

guppy_basecaller -i $INPUT -s $FOLDER -c dna_r9.4.1_450bps_sup.cfg --device auto --recursive --detect_mid_strand_adapter --min_qscore 7 --barcode_kits SQK-NBD110-24
```
```cat Guppy_Super_Accurate_model_23EV612/*.fastq > 23EV612_Guppy.fastq ```

#### Filtlong
[Filtlong](https://github.com/rrwick/Filtlong) was used to remove reads shorter than 1kbp and exclude the worst 5% of the reads. 
```bash
module load Filtlong/0.2.0
filtlong --min_length 1000 --keep_percent 95 23EV612_Guppy.fastq| gzip > 23EV612_Guppy_Filtlong.fastq.gz
```
#### Porechop
[Porechop](https://github.com/rrwick/Porechop) was used to remove any adaptors 
```bash
#!/bin/bash -e
#SBATCH -c 15 --mem=16Gb --time 00:120:00 -J Porechop_EV

module purge
module load Porechop/0.2.4-gimkl-2020a-Python-3.8.2

for input in *_Guppy_Filtlong.fastq.gz
do
base=$(basename ${input} _Guppy_Filtlong.fastq.gz)
echo "working with file $input"
echo "basename is $base"

OUTPUT=/pathtoprocessedreads/StaphA_${base}_Guppy_Filtlong_Porechop.fastq.gz
porechop -i $input -o $OUTPUT --discard_middle

done
```
#### NanoStat
[NanoStat](https://github.com/wdecoster/nanostat) was used to perform quality checks to determine the mean read quality
```
module load NanoStat/1.5.0-gimkl-2020a-Python-3.8.2
NanoStat --fastq StaphA_23EV612_Guppy_Filtlong_Porechop.fastq.gz
```
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
```
module load chopper/0.2.0-GCC-11.3.0
gunzip -c StaphA_23EV612_Guppy_Filtlong_Porechop.fastq.gz | chopper -q 18 | gzip > StaphA_23EV612_Guppy_Filtlong_Porechop_Chopper.fastq.gz
#The q is the mean read quality as specified from NanoStat.
```
#### Flye 
[Flye](https://github.com/mikolmogorov/Flye) was used as the *de novo* assembler to assemble the isolate 23EV612. 
```bash
#!/bin/bash -e
#SBATCH --nodes 1 --cpus-per-task 1 --ntasks 10 --mem=50G -J Flye_EV --time=72:00:00 --hint=nomultithread

module load Flye/2.9.1-gimkl-2022a-Python-3.10.5
INPUT=StaphA_23EV612_Guppy_Filtlong_Porechop_Chopper.fastq.gz
OUTPUT=23EV612_Flye/
flye --nano-hq $INPUT --genome-size 2.8m -o $OUTPUT -t 10 -i 3 --asm-coverage 50
```
#### Medaka 
[Medaka](https://github.com/nanoporetech/medaka) was used as the polishing tool for the 23EV612. 
```bash
#!/bin/bash -e
#SBATCH --nodes 1 --cpus-per-task 1 --job-name medaka_Flye --mem=20G --time=24:00:00 --output=%x_%j.out --error=%x_%j.err --hint=nomultithread
module load medaka/1.6.0-Miniconda3-4.12.0

READS=StaphA_23EV612_Guppy_Filtlong_Porechop_Chopper.fastq.gz
CONTIG_FILE=23EV612_Flye/assembly.fasta

medaka_consensus -i $READS -d $CONTIG_FILE -o medaka_23EV612/ -t 2 -m r941_min_sup_g507
```

#### Bwa mem and samtools sort 

#### Pilon

#### CheckM Quast 












## Nullarbor
[Nullarbor](https://github.com/tseemann/nullarbor) is the genomic analysis pipeline chosen for bovine, human and ST1 *S. aureus* analysis. In the nullarbor.pl file --mincov was specified as 60 and --minid as 90. The script was run on NESI. 

```bash
#!/bin/bash -e
#SBATCH --cpus-per-task=20 --mem 160Gb --time 166:00:00 -J nullarbor_EV
#Alter Reference Genome for human, bovine or ST1 analysis - comment out the one not needed and alter Reference in the script.
#Alter text file for purpose - bovine, human or ST1 analysis
Ref_Bovine = SaureusRF122.fa
#Ref_Human = SaureusMSSA476.fa
#Ref_ST1 = 23EV612.fa

mkdir -p $TMPDIR
module purge
module load nullarbor/2.0.20191013

nullarbor.pl --name StaphAureusIsolates --ref $Ref_Bovine \
--input BovineIsolates.txt \
--outdir BovineAnalysis --force \
--cpus $SLURM_CPUS_PER_TASK --run --mlst saureus --trim --taxoner-opt confidence 0.5 --treebuilder iqtree --treebuilder-opt "-wbtl -st DNA -b 1000 --alrt 1000" --annotator-opt "--fullanno --compliant"
```

## Defining *Spa* Types 
*Spa* Types were determined using [SpaTyper](https://github.com/HCGB-IGTP/spaTyper). The list of genomes is made from the contigs.fa created by Nullarbor. 
```bash 
ls *.fa > list_of_genomes.txt ##Creating a list of genomes from all fastq files 
sed -i '1 i\spaTyper  -f ' list_of_genomes.txt
echo "--output spa.xls" &>> list_of_genomes.txt
tr '\n' ' ' < list_of_genomes.txt > spa_command.txt
chmod +x spa_command.txt
./ spa_command.txt
```

## Defining *agr* types
[AGRvate](https://github.com/VishnuRaghuram94/AgrVATE) was used to identify each genomes *agr* type. A conda environment was created then the contig fasta files was used as input for AGRvate. 

## Hybrid Genome

## Gubbins 


## IQtree




## Python script for SNP distance matrix 

## AMR_stats.py 

## Virulence Gene Analysis 
Virulence Genes were split into their broad function: Adherence, Enterotoxin, Exoenzyme, Exotoxin, Haemolysin, Immune Modulation, Intracellular Adhesion, Type VII Secretion System and Others. The csv file was a binary matrix where Isolate in column 1, host in column 2 and the following columns with the genes detected where 1 is present and 0 being absent. Each individual function csv file was the input into the below python script to put *S. aureus* isolates of each host  into groups to identify the combination of genes that are present in majority of each bovine and human *S. aureus* host. 

```python
file_path = "exoenzyme.csv" #Altered for each functional group
results = {}

with open(file_path,'r') as data:
    for line in data.readlines():
        if 'Source' in line:
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
    ##To run python scriptname > outputfile
```

