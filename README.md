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

mkdir -p $TMPDIR
module purge
module load nullarbor/2.0.20191013

nullarbor.pl --name StaphAureusIsolates --ref $Ref_Bovine \
--input BovineIsolates.txt \
--outdir BovineAnalysis --force \
--cpus $SLURM_CPUS_PER_TASK --run --mlst saureus --trim --taxoner-opt confidence 0.5 --treebuilder iqtree --treebuilder-opt "-wbtl -st DNA -b 1000 --alrt 1000" --annotator-opt "--fullanno --compliant"
```
## FastQC 

## SpaTyper 

## AGRvate 

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

