Guppy
```
#!/bin/bash -e
#SBATCH -J guppy_gpu --gpus-per-node A100:1 --mem 6G --cpus-per-task 4 --time 10:00:00 --output slurmout.%j.out

module purge
module load ont-guppy-gpu/6.4.6

INPUT=Path/to/raw/reads/20191112_0850_MN30056_FAN36754_626723d4/
FOLDER=Guppy_Super_Accurate_model_23EV612/

guppy_basecaller -i $INPUT -s $FOLDER -c dna_r9.4.1_450bps_sup.cfg --device auto --recursive --detect_mid_strand_adapter 
```
```
cat $FOLDER/*.fastq > 23EV612.fastq 
```

Filtlong 
```
module load Filtlong/0.2.0
filtlong --min_length 1000 --keep_percent 95 23EV612_Guppy.fastq| gzip > 23EV612_Guppy_Filtlong.fastq.gz
```

Porechop
```
#!/bin/bash -e
#SBATCH -c 15 --mem=16Gb --time 00:120:00 -J Porechop_EV

module purge
module load Porechop/0.2.4-gimkl-2020a-Python-3.8.2

for input in *_Guppy_Filtlong.fastq.gz
do
base=$(basename ${input} _Guppy_Filtlong.fastq.gz)
echo "working with file $input"
echo "basename is $base"

OUTPUT=StaphA_${base}_Guppy_Filtlong_Porechop.fastq.gz
porechop -i $input -o $OUTPUT --discard_middle

done
```
NanoStat 
```
module load NanoStat/1.5.0-gimkl-2020a-Python-3.8.2
NanoStat --fastq StaphA_23EV612_Guppy_Filtlong_Porechop.fastq.gz
```
Chopper
```
module load chopper/0.2.0-GCC-11.3.0
gunzip -c StaphA_23EV612_Guppy_Filtlong_Porechop.fastq.gz | chopper -q 18 | gzip > StaphA_23EV612_Guppy_Filtlong_Porechop_Chopper.fastq.gz
#The q is the mean read quality as specified from NanoStat.
```
Flye
```
#!/bin/bash -e
#SBATCH --nodes 1 --cpus-per-task 1 --ntasks 10 --mem=50G -J Flye_EV --time=72:00:00 --hint=nomultithread

module load Flye/2.9.1-gimkl-2022a-Python-3.10.5
INPUT=StaphA_23EV612_Guppy_Filtlong_Porechop_Chopper.fastq.gz
OUTPUT=23EV612_Flye/
flye --nano-hq $INPUT --genome-size 2.8m -o $OUTPUT -t 10 -i 3 --asm-coverage 50
```
Medaka
```
#!/bin/bash -e
#SBATCH --nodes 1 --cpus-per-task 1 --job-name medaka_Flye --mem=20G --time=24:00:00 --output=%x_%j.out --error=%x_%j.err --hint=nomultithread
module load medaka/1.6.0-Miniconda3-4.12.0

READS=StaphA_23EV612_Guppy_Filtlong_Porechop_Chopper.fastq.gz
CONTIG_FILE=23EV612_Flye/assembly.fasta

medaka_consensus -i $READS -d $CONTIG_FILE -o medaka_23EV612/ -t 2 -m r941_min_sup_g507
```
```
#Rename medaka polished output
mv medakaoutput.fasta 23EV612_Medaka_Polish.fasta
```
BWA & SAMtools
```
#!/bin/bash -e
#SBATCH -c 10 --mem=8Gb --time 00:10:00 -J bwamem_EV

module load BWA/0.7.17-GCC-11.3.0
module load SAMtools/1.16.1-GCC-11.3.0

REF=23EV612_Medaka_Polish.fasta
bwa mem -t 6 -R"@RG\tID:23EV612\tPL:ILLUMINA\tSM:LIC_23EV612_2019" $REF 23EV612_R1.fq.gz 23EV612_R2.fq.gz | samtools view - -Sb | samtools sort - -@10 -o 23EV612_Medaka_Illumina.bam
```
```
#The output bam is indexed using SAMtools 
module load SAMtools/1.16.1-GCC-11.3.0
samtools index 23EV612_Medaka_Illumina.bam
```
Pilon & BWAmem_Pilon
```
#!/bin/bash -e
#SBATCH -c 4 --mem=10Gb --time 00:10:00 -J PILON_EV

mkdir -p 23EV612_Pilon/
OUTPUT=23EV612_Pilon/ 

module load Pilon/1.24-Java-15.0.2

GENOME=23EV612_Medaka_Polish.fasta
BAM=23EV612_Medaka_Illumina.bam 

java -Xmx10G -jar $EBROOTPILON/pilon.jar --genome $GENOME --fix all --changes --frags $BAM --threads 10 --output $OUTPUT/pilon_round1 | tee $OUTPUT/round1.pilon
```
```
#The output fasta file is indexed using BWA
module load BWA/0.7.17-GCC-11.3.0
samtools index pilon_round1.fasta 
```
```
#!/bin/bash -e
#SBATCH -c 5 --mem=8Gb --time 00:10:00 -J bwamempilon_EV

module load BWA/0.7.17-GCC-9.2.0
module load SAMtools/1.16.1-GCC-11.3.0

INPUT=pilon_round1.fasta
ILLUMINAR1=23EV612_R1.fq.gz
ILLUMINAR2=23EV612_R2.fq.gz
bwa mem -t 5 $INPUT $ILLUMINAR1 $ILLUMINAR2 | samtools view - -Sb | samtools sort - -@14 -o Pilon1_23EV612.bam
```bash
The output bam is indexed using SAMtools 
module load SAMtools/1.16.1-GCC-11.3.0
samtools index Pilon1_23EV612.bam
```
CheckM & QUAST
```
#!/bin/bash -e
#SBATCH -c 5 --mem=10Gb --time 00:30:00 -J CHECKM_EV
module load CheckM/1.2.1-gimkl-2022a-Python-3.10.5
module load QUAST/5.2.0-gimkl-2022a

OUTPUT=CheckM_Quast/
INPUT=PolishedGenomes/

checkm taxon_set genus "Staphylococcus" staph.ms 
checkm analyze -t 16 -x fa staph.ms $INPUT $OUTPUT/checkmoutput 
checkm qa -t 16 staph.ms $OUTPUT/checkmoutput
##Need to review contamination and completeness - contamination less than 5% and completeness above 90%

quast.py $INPUT/*.fa
```
