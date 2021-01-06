# CREOLE
Reference-free identification of open reading frames and encoded proteins from nanopore transcriptomic long reads. 
# Requirements
- Minimap2
- GeneMark
- DIAMOND
# Installation
- Clone the repository

`git clone --recurse-submodules https://github.com/comprna/RATTLE`

# Usage example 
**Warning** All the commands and parameters are still highly experimental and subject to changes in future versions

- in cluster (Slurm)

```
#!/bin/bash
#SBATCH --partition=normal
#SBATCH --nodes=4
#SBATCH -J creoleT4_RNA_mecat
#SBATCH --mem=91G
#SBATCH -c 24
#SBATCH -e error/creoleT4_RNA_mecat.err
#SBATCH -o error/creoleT4_RNA_mecat.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=joelamaque@gmail.com

module load Python

export PATH=/users/genomics/joel/creoleT4:$PATH
export PATH="$PATH:/genomics/users/joel/soft/ST_alg/gmst_linux_64"
export PATH="$PATH:/projects_rg/joel/tools/diamond"

dir="/creole/RNA/mecat"
mkdir -p $dir
cd $dir

array=(`ls /projects_rg/joel/charlote_et_all/nanopore/corr_reads/RNA/mecat/*_cut_150_cor.fa_mix | cut -d"/" -f9 | cut -d"." -f1`)
array1=(`ls /projects_rg/joel/charlote_et_all/nanopore/corr_reads/RNA/mecat/*_cut_150_cor.fa_mix`)

for index in ${!array[@]}; do
        echo ${array[$index]} ${array1[$index]}
        python /users/genomics/joel/creoleT3/creole.py\
        -gmst /genomics/users/joel/soft/ST_alg/gmst_linux_64/gmst.pl\ 
        -reads ${array1[$index]}\ 
        -o ${array[$index]}\  
        -ct prok\ 
        -st direct\ 
        -diamond /projects_rg/joel/tools/diamond\ 
        -diamond_bd /users/genomics/joel/charlote_et_all/orfs/prot_db/nr.dmnd\ 
        -f 2\ 
        -minimap2 /projects_rg/joel/tools/minimap2-2.17_x64-linux/minimap2\ 
        -mr /projects_rg/joel/gencode_v31/gencode.v31.transcripts_noPseudogene.fa\ 
        -prot /projects_rg/joel/gencode_v31/gencode.v31.pc_translations.fa
done
```

- Sample commands
```
python creole.py -gmst gmst.pl -reads read.fa -o file.out -ct prok -st direct -diamond diamond -diamond_bd nr.dmnd -f 2 -minimap2 minimap2 -mr gencode.v31.transcripts.fa -prot gencode_v31/gencode.v31.pc_translations.fa

```
- Commands and options

**-reads:** input fasta/fastq file

**-gmst:** Path to the GeneMark

**-diamond:** Path to the DIAMOND

**-diamond_bd:** protain database

**-minimap2:** Path to the Minimap2

**-mr:** transcript reference

**-prot:** protain reference

**-o:** output files

**-ct:** cell type (prok or euka)

**-st:** strand (direct, reverse and both)

**-f:** types of filters to select the best read match to transcript (1 - my score, 2 - mapping quality score, 3 - max mapping quality (60) and 4 - DP alignment score)
