# Comparative Genomics ST1 Staphylococcus aureus

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
nullarbor.pl --name StaphAureusIsolates --ref $Ref_Bovine \
--input BovineIsolates.txt \
--outdir BovineAnalysis --force \
--cpus $SLURM_CPUS_PER_TASK --run --mlst saureus --trim --taxoner-opt confidence 0.5 --treebuilder iqtree --treebuilder-opt "-wbtl -st DNA -b 1000 --alrt 1000" --annotator-opt "--fullanno --compliant"
```
## VFG Headers

## FastQC 

## SpaTyper 

## AGRvate 

## Hybrid Genome

## Gubbins 

## IQtree

## Python script for SNP distance matrix 

## AMR_stats.py 
